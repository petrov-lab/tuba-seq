#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import os
import argparse
import params
from tuba_seq import fastq
from tuba_seq.shared import logPrint

############### Input Parameters that will be retained ########################

cluster_flank = 7           # Parameters used in Rogers et al (2017) Nat meth.
training_flank = 17         # These have changed in the latest version. 

############### Input Parameters that will be deprecated ######################

KMERS = 2                   # Number of kmer searches used to find the beginning and end of double-barcodes
                            # Hard-coded as 2, do not change.
match_len = 6               # Length of each kmer
symmetric_immediate_truncation_of_read = slice(len(params.head) - len(params.tail), len(params.master_read)) 
    # Unused sections of reads are immediately truncated to accelerate processing. 
    # This truncation assumes len(head) > len(tail)

###############################################################################

truncated_master_read = params.master_read[symmetric_immediate_truncation_of_read]
cluster_distance_from_start = params.barcode_length + cluster_flank

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.")
parser.add_argument("--base_dir", type=str, default=os.getcwd(),
help='Base directory to work within. This directory must contain a folder entitled "{:}" containing all FASTQ files to process.'.format(params.original_dir))
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multithreaded operation')
parser.add_argument('-e', '--maxEE', type=float, default=2, help='Maximum Expected Errors per read (default 2.)')

args = parser.parse_args()
base_dir = args.base_dir.rstrip('/')
os.chdir(base_dir)

Log = logPrint(verbose=args.verbose)        # An object to save all of the program's output based on the choice of verbosity

if args.verbose and args.parallel:
    print("Verbose output is incompatible with parallel operation. Will use single thread...")
    args.parallel = False

if args.parallel:
    from tuba_seq.pmap import pmap as map

for Dir in [params.preprocessed_dir, params.training_dir]:
    if not os.path.exists(Dir):
        os.makedirs(Dir)

files = list(filter(lambda fn: params.fastq_handle in fn, os.listdir(params.original_dir)))

Log('Processing {:} Files in {:}/original/.'.format(len(files), base_dir))

def easy_N_fixes(DNA, searcher=re.compile('(...N...)'), maxN=3):
    """Repair N bases in DNA using params.master_read. 

Keyword Arguments:
searcher -- the regular expression used to search for N bases
maxN -- the maximum number of N bases to try to repair (default 3)

This function cannot always repair N bases, e.g. in barcode regions.
    """
    
    if DNA.count('N') > maxN:
        return DNA
    broken = searcher.split(DNA)
    try:
        return ''.join([re.findall(s.replace('N', '.'), truncated_master_read)[0] if 'N' in s else s for s in broken])
    except IndexError:
        return DNA

offsets = match_len*np.arange(KMERS)
head_matchers = [fastq.singleMismatcher(params.head[-i-match_len:][:match_len]) for i in offsets]
tail_matchers = [fastq.singleMismatcher(params.tail[i:i+match_len]) for i in offsets]
    
start_expected = min(len(params.head), len(params.tail))
stop_expected = start_expected + params.barcode_length

def process_file(f):
    """Processes a FASTQ file into a training file for DADA2 & a barcode clustering file
"""
    df = (fastq.fastqDF.from_file(os.path.join(params.original_dir, f), use_Illumina_filter=True, fake_header=False)
                       .co_slice(symmetric_immediate_truncation_of_read))
    short_filename = f.split(params.fastq_handle)[0]
    
    # Attempt #1 to repair N bases
    degen = df['DNA'].str.contains("N")
    problems = df.loc[degen, 'DNA']
    df.loc[degen, 'DNA'] = problems.apply(easy_N_fixes)
    
    # Identify the beginning of the barcode region 
    heads = [df['DNA'].apply(func.find) for func in head_matchers]
    starts = heads[0] + match_len
    wrong_starts = (starts - start_expected).abs() > params.allowable_deviation
    initial_wrongs = wrong_starts.sum()
    starts[wrong_starts] = heads[1][wrong_starts] + match_len*2
    wrong_starts = (starts - start_expected).abs() > params.allowable_deviation
    final_wrongs = wrong_starts.sum()
    
    # Identify the end of the barcode region
    tails = [df['DNA'].apply(func.find) for func in tail_matchers]
    stops = tails[0]
    wrong_stops = (stops - stop_expected).abs() > params.allowable_deviation
    truncated_degen = ((stops > starts) & (stop_expected - stops > params.allowable_deviation)).sum()
    Log("{:.1%} of reads had truncated/missing sgIDs & barcodes.".format(truncated_degen/len(df)))
    initial_wrongs += wrong_stops.sum() 
    stops[wrong_stops] = tails[1][wrong_stops] - match_len
    wrong_stops = (stops - stop_expected).abs() > params.allowable_deviation
    final_wrongs += wrong_stops.sum()
   
    # Discard reads with start/stop of barcode in the wrong location
    keep = ~(wrong_starts | wrong_stops)
    passed_kmers = df.select_reads(keep)
    Log('{:.1%} passed kmer tests.'.format(len(passed_kmers)/len(df)))
    starts = starts[keep]
    stops = stops[keep]

    # Only train DADA2 on reads lacking N bases, even if we can repair these N bases.
    # Otherwise, the error rate estimate will be biased. 
    was_not_degen = ~passed_kmers['QC'].str.contains('#')
    train = passed_kmers.select_reads(was_not_degen)
    Log('Error training will use {:.1%} of remaining reads.'.format(len(train)/len(passed_kmers)))
   
    # Train DADA2 only on the barcode-flanking regions of reads
    train_starts = starts[was_not_degen].apply(lambda start: slice(start-training_flank, start))
    train_stops  =  stops[was_not_degen].apply(lambda stop:  slice(stop, stop + training_flank))
    (train.vector_slice(train_starts, train_stops)
          .drop_abnormal_lengths()
          .write(params.training_dir+f))

    # Trim reads for DADA2 clustering 
    slices = starts.apply(lambda start: slice(start - cluster_flank, start + cluster_distance_from_start))
    sliced = passed_kmers.vector_slice(slices).drop_abnormal_lengths()
    cluster = (sliced.select_reads(~sliced.isDegenerate())
                     .drop_abnormal_lengths())
    
    # Discard reads with poor QC scores
    EEs = cluster.expected_errors()
    clean = cluster.select_reads(EEs <= args.maxEE)

    # Write clustering file
    clean.write(params.preprocessed_dir+('{short_filename}.fastq'.format(**locals())))
    
    #Output summary
    counts = dict(saved = initial_wrongs - final_wrongs,
                  reads = len(df),
           passed_kmers = len(passed_kmers),
                cluster = len(cluster),
                  clean = len(clean),
              truncated = truncated_degen,
            good_flanks = len(passed_kmers)+truncated_degen,
                    EEs = EEs.sum())

    reads = counts['reads']
    percents = {k:v/reads for k, v in counts.items()}
    percents['EE_rate'] = counts['EEs']/(len(clean)*len(clean['DNA'].values[0]))
    Log("{:} ({:.2}M reads): {good_flanks:.0%} good flanks, {clean:.0%} used. Estimated Error Rate: {EE_rate:.3%}.".format(
        short_filename, reads*1e-6, **percents), print_line=True)
    return counts

all_output = map(process_file, files)
output_df = pd.DataFrame(list(all_output), index=files)

totals = output_df.sum()

reads = totals['reads']
mean_error_rate = totals.pop('EEs')/totals['clean']/(params.barcode_length + 2*cluster_flank)

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
            maxEE=args.maxEE,
            error_rate=mean_error_rate,
            **percents), print_line=True)
Log.close()

