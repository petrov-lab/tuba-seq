import pandas as pd
import numpy as np
import re
import os
import argparse
from fastq import fastqDF, singleMismatcher 
from params import master_read, head, tail, core_length


############### Input Parameters that will be retained ########################

from params import cluster_flank, training_flank, maxEE, allowable_deviation, preprocessed, training, fastq_handle
#cluster_flank = 7           # Parameters used in Rogers et al; subject to change
#training_flank = 17

############### Input Parameters that will be deprecated ######################

kmers = 2                   # Number of kmer searches used to find the beginning and end of double-barcodes
match_len = 6               # Length of each kmer
symmetric_immediate_truncation_of_read = slice(len(head) - len(tail), len(master_read)) # Unused sections of reads are immediately truncated to accelerate processing. 
                                                                                        # This truncation assumes len(head) > len(tail)
from logPrint import logPrint                                                           # Currently, output goes to stdout and a log file; this will change
Log = logPrint('preprocessing_Log.txt')                                                 
###############################################################################

truncated_master_read = master_read[symmetric_immediate_truncation_of_read]
cluster_distance_from_start = core_length + cluster_flank

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.")
parser.add_argument("--base_dir", type=str, help='Base directory to work within. This directory must contain a folder entitled "original" containing all FASTQ files to process.', default=os.getcwd())
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Parallelize operation')
args = parser.parse_args()
base_dir = args.base_dir.rstrip('/')
os.chdir(base_dir)

try: 
    import pmap
except ImportError:
    if args.parallel:
        print("Cannot import 'multiprocessing' module. Parallelization not possible.")
else:
    if args.parallel:
        map = pmap.pmap

sg_info = pd.read_csv(args.sgRNA_file)
for Dir in [preprocessed, training]:
    if not os.path.exists(Dir):
        os.makedirs(Dir)

files = list(filter(lambda fn: fastq_handle in fn, os.listdir('original/')))

Log('Processing {:} Files in {:}/original/.'.format(len(files), base_dir))

def easy_N_fixes(DNA, searcher=re.compile('(...N...)'), maxN=3):
    if DNA.count('N') > maxN:
        return DNA
    broken = searcher.split(DNA)
    try:
        return ''.join([re.findall(s.replace('N', '.'), truncated_master_read)[0] if 'N' in s else s for s in broken])
    except IndexError:
        return DNA

offsets = match_len*np.arange(kmers)
head_matchers = [singleMismatcher(head[-i-match_len:][:match_len]) for i in offsets]
tail_matchers = [singleMismatcher(tail[i:i+match_len]) for i in offsets]
    
start_expected = min(len(head), len(tail))
stop_expected = start_expected + core_length

def process_file(f):
    short_f = f.split(fastq_handle)[0]
    df = fastqDF('original/'+f, symmetric_immediate_truncation_of_read)
    degen = df.isDegenerate()
    problems = df.loc[degen, 'DNA']
    df.loc[degen, 'DNA'] = problems.apply(easy_N_fixes)

    heads = [df['DNA'].apply(func.find) for func in head_matchers]
    tails = [df['DNA'].apply(func.find) for func in tail_matchers]
   
    starts = heads[0] + match_len
    wrong_starts = (starts - start_expected).abs() > allowable_deviation
    initial_wrongs = wrong_starts.sum()
    starts[wrong_starts] = heads[1][wrong_starts] + match_len*2
    wrong_starts = (starts - start_expected).abs() > allowable_deviation
    final_wrongs = wrong_starts.sum()
    stops = tails[0]
    wrong_stops = (stops - stop_expected).abs() > allowable_deviation
    truncated_degen = ((stops > starts) & (stop_expected - stops > allowable_deviation)).sum()
    
    if args.verbose:
        Log("{:.1%} of reads had truncated/missing sgIDs & barcodes.".format(truncated_degen/len(df)))
    
    initial_wrongs += wrong_stops.sum() 
    stops[wrong_stops] = tails[1][wrong_stops] - match_len
    wrong_stops = (stops - stop_expected).abs() > allowable_deviation
    final_wrongs += wrong_stops.sum()
    keep = -(wrong_starts | wrong_stops)
    passed_kmers = df.loc[keep, :]
    passed_kmers.__class__ = fastqDF

    if args.verbose:
        Log('{:.1%} passed kmer tests.'.format(len(passed_kmers)/len(df)))
    
    starts = starts[keep]
    stops = stops[keep]

    was_not_degen = -passed_kmers['QC'].str.contains('#')
    train = passed_kmers.loc[was_not_degen, :]
    train.__class__ = fastqDF
    
    if args.verbose:
        Log('Error training will use {:.1%} of remaining reads.'.format(len(train)/len(passed_kmers)))
    
    train_starts = starts[was_not_degen].apply(lambda start: slice(start-training_flank, start))
    train_stops  =  stops[was_not_degen].apply(lambda stop:  slice(stop, stop + training_flank))

    train = train.vslice(train_starts, train_stops, enforce=True) 
    train.write(training+f)

    slices = starts.apply(lambda start: slice(start - cluster_flank, start + cluster_distance_from_start))
    sliced = passed_kmers.vslice(slices, enforce=True)

    cluster = sliced.loc[-sliced.isDegenerate(), :]
    cluster.__class__ = fastqDF
    cluster.drop_abnormal_lengths()
    
    try:
        clean, EEs = cluster.expected_errors(maxEE=maxEE)
    except ValueError as e:
        print(short_f, 'had an invalid Q score.')
        raise e

        #OUTPUT
    clean.write(preprocessed+('{short_f}.fastq'.format(**locals())))
    
    counts = dict(saved = initial_wrongs - final_wrongs,
                  reads = len(df),
           passed_kmers = len(passed_kmers),
                cluster = len(cluster),
                  clean = len(clean),
              truncated = truncated_degen,
            good_flanks = len(passed_kmers)+truncated_degen,
                    EEs = EEs)

    reads = counts['reads']
    percents = {k:v/reads for k, v in counts.items()}
    Log("{:} ({:.2}M reads): {good_flanks:.0%} good flanks, {clean:.0%} used, {unknown_IDs:.0%} unknown sgIDs. Estimated Error Rate: {:.3%}.".format(short_f, reads*1e-6, EEs/len(clean)/(core_length + 2*cluster_flank), **percents))
    return counts

all_output = map(process_file, files)
output_df = pd.DataFrame(all_output)

totals = output_df.sum()

reads = totals['reads']
mean_error_rate = totals.pop('EEs')/totals['clean']/(core_length + 2*cluster_flank)

percents = {k:v/reads for k, v in totals.items()}

Log("""
Summary for {:}:
---------------------------------------
Processed {:.2}M Reads.
{good_flanks:.1%} had reasonable headers & tails ({saved:.1%} were saved by 2nd kmer).
{passed_kmers:.1%} had sgIDs & barcodes.
{cluster:.1%} had fixable Ns.
{clean:.1%} had <= {maxEE:g} expected errors and the average Error Rate was {error_rate:.4%}.
""".format( base_dir, 
            reads*1e-6, 
            maxEE=maxEE,
            error_rate=mean_error_rate,
            **percents))
Log.close()

