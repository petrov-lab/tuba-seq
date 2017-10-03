#!/usr/bin/env python3
import pandas as pd
import os, numpy, argparse
from tuba_seq.fastq import fastqDF, infer_master_read, MasterRead
from tuba_seq.shared import logPrint
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from collections import defaultdict
import warnings

fastq_ext = '.fastq'
############### Input Parameters ########################

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input_dir', type=str, 
    help='Input directory with fastq files.') # By default will be {:} if using merged paired-end reads and {:} if not.'.format(merged_dir, original_dir))
parser.add_argument('--master_read', type=str, default='infer', help="Outline of the amplicon sequence, degenerate bases can be either 'N' or '.' E.g. ATCGATCGTCGA........AA.....TT.....AA.....GTTCGAGCTAGC")

parser.add_argument('-t', '--training_dir', default='training/', help='Directory to save files for error training.') 
parser.add_argument('-o', '--output_dir', default='preprocessed/', help='Directory to save files for barcode clustering.') 
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multi-threaded operation')
parser.add_argument('-s', '--search_blast', action='store_true', help='Use NCBI BLAST algorithm to identify contaminations in samples')
parser.add_argument('-l', '--local_blast', action='store_true', help='Use local NCBI BLAST+ algorithm to accelerate searching (if present), see tuba_seq/blast.py.')
parser.add_argument('-d', '--derep', action='store_true', help='De-replicate output fastQ files for DADA2 to minimize file sizes.')
parser.add_argument('-f', '--fraction', type=float, default=0.05, help='Minimum fraction of total reads to elicit a BLAST-search of an unknown sequence.')
parser.add_argument('-k',  '--skip', action='store_true', help='Skip files that already exist in output directories.')

parser.add_argument('-a', '--allowable_deviation', type=int, default=4, help="Length of Indel to tolerate before discarding reads.")
parser.add_argument('--alignment_flank', type=int, default=22, help='# of bases flanking the degenerate region to be used to score the quality of the read.')
parser.add_argument('--training_flank', type=int, default=18, help='# of bases flanking the degenerate region to be used to develop the DADA2 error model.')
parser.add_argument('-M', '--min_align_score', type=float, default=0.6, help='Minimum alignment score needed to keep read, Range [0, 1).')
parser.add_argument('--compression', default='bz2', choices=['bz2', 'gzip', 'none'], help='Compression algorithm for saved file.')
###############################################################################
args = parser.parse_args()

output_dirs = args.training_dir, args.output_dir
dada2 = importr("dada2")
R_base = importr('base')

if args.derep or args.compression == 'none':
    compression_ext = ''
elif args.compression == 'gzip':
    from gzip import open       
    compression_ext = '.gz'
elif args.compression == 'bz2':
    from bz2 import open
    compression_ext = '.bz2'

if args.parallel:
    from tuba_seq.pmap import low_memory_pmap as map 

for Dir in output_dirs:
    os.makedirs(Dir, exist_ok=True)

files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if fastq_ext in f]
max_file_size = 1e9 if args.parallel else 2e9


Log = logPrint(args)
Log('Processing {:} samples found in {:}.'.format(len(files), args.input_dir), print_line=True)
def process_sample(filename):
    sample = os.path.basename(filename.partition(fastq_ext)[0])
    filenames = [os.path.join(Dir, sample)+fastq_ext+compression_ext for Dir in output_dirs]
    if args.skip and all(map(os.path.isfile, filenames)):
        Log("Skipping "+sample+" (already exists).")
        return None
    if os.path.getsize(filename) > max_file_size:
         warnings.warn("{:} is {:.1} GB. preprocess.py must load this entire fastq file into memory--and gzip decompression will inflate this size ~10x. You'll need about twice this amount of memory for successful execution.", RuntimeWarning)
         if args.parallel:
            warnings.warn("Running preprocess.py in non-multi-threaded form saves memory (multi-threading processes several fastq files simultaneously).", RuntimeWarning)

    reads = fastqDF.from_file(filename, fake_header=True, use_Illumina_filter=True)
    init_reads = len(reads)
    gb = reads.set_index("DNA")['QC'].groupby(level='DNA') 
    counts = gb.count()
    

    master_read = MasterRead(args.master_read if args.master_read != 'infer' else infer_master_read(counts.nlargest(1000).index), 
                             args.allowable_deviation, args.min_align_score, args.alignment_flank, args.training_flank, args.allowable_deviation)
    with open(filenames[0], 'wb') as training_file, open(filenames[1], 'wb') as cluster_file:
        reads_tuples = gb.agg(master_read.process_read_set, reads.construct_read_set, training_file, cluster_file)
    output = pd.DataFrame(reads_tuples.tolist(), columns=['outcome', 'score'], index=counts.index)
    output['Reads'] = counts.values
    unaligned = output.query('outcome == "Unaligned"')['Reads']
    PhiX = pandas2ri.ri2py(dada2.isPhiX(pandas2ri.py2ri(unaligned.index))) == 1
    contaminants = unaligned.loc[(~PhiX)&(unaligned>int(args.fraction*init_reads))]
    tallies = output.groupby('outcome')['Reads'].sum()
    tallies['PhiX'] = unaligned.loc[PhiX].sum()
    tallies['Unaligned'] = tallies['Unaligned'] - tallies['PhiX'] if 'Unaligned' in tallies else 0
    percents = defaultdict(float, tallies/init_reads)
    Log('Sample {:}: {:.2f}M Reads, {Unaligned:.1%} unaligned, {PhiX:.1%} PhiX, {Wrong Barcode Length:.1%} inappropriate barcode length, {Clustered:.1%} Cluster-able.'.format(sample, init_reads*1e-6, **percents))
    tallies['Initial'] = init_reads
    reads.info['Master Read'] = master_read.ref
    return dict(tallies=tallies,
                info=reads.info, 
                unknown_DNAs=contaminants/init_reads,
                scores=output.groupby('score').sum())

_reports = list(map(process_sample, files))

dfs = {r['info'].loc['Sample']:r for r in _reports if r is not None}
meta_datas = ['tallies', 'unknown_DNAs', 'info', 'scores']
output = {meta_data:pd.concat({d['info'].loc['Sample']:d[meta_data] for d in _reports if d is not None}, names=['Sample']) for meta_data in meta_datas}

def derep(filename):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        try:
            derep = dada2.derepFastq(filename, verbose=args.verbose)
            R_base.saveRDS(derep, file=filename.replace(fastq_ext, '.rds'))
        except Exception: 
            print("Could not derep "+filename)
        else:
            os.remove(filename)

if args.derep:
    Log("De-replicating outputs for DADA2...")
    list(map(derep, [os.path.join(Dir, filename) for Dir in output_dirs for filename in os.listdir(Dir) if fastq_ext in filename]))

if args.search_blast and len(output['unknown_DNAs']) > 0: 
    from tuba_seq.blast import sleuth_DNAs
    gb = output['unknown_DNAs'].groupby(level='DNA')
    DNAs = list(gb.keys())
    Log("BLAST-searching", len(DNAs), "common unknown sequence(s)...", print_line=True)
    searches = sleuth_DNAs(DNAs, local_blast=args.local_blast)
    if len(searches) == 0:
        Log("Could not find any acceptable matches.", print_line=True)
    else:
        Log("                -- Best Match --                          | E-score | Sample(s) | % of Sample", print_line=True)
        for dna, row in searches.iterrows():
            group = gb.groups[dna]
            blast_result = "{short_title:<60} {e:1.0e}".format(**row)
            for (_dna, sample), fraction in group.iteritems():
                Log("{:}    {:}    {0:.2%}".format(blast_result, sample, fraction), print_line=True)
                blast_result = len(blast_result)*' '

tallies = output['tallies'].astype(int).unstack()
info = output['info'].unstack()

if args.master_read == 'infer':
    MRs = info['Master Read']
    if len(MRs.value_counts()) == 1:
        Log("Inferred the Master Read to be:\n"+MRs[0], True)
    else:
        Log("""Did not infer the same Master Read from every Sample!!
Inferred the following Master Reads:
{:}
**You probably want to specify --master_read explicitly**""".format(MRs.to_string()), True)

if len(info['Instrument'].value_counts()) > 1:
    Log("This run contains fastq files from two different Illumina machines. This is not recommended.", print_line=True)

def merge_with_metadata_file(summary):
    pass

def plot_alignment_histogram(data, filename='alignment_histogram.pdf', hist_kwargs={'bins':numpy.linspace(0, 1, 40)}):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    ax.hist(data.index.values, weights=data.values, **hist_kwargs)
    ax.axvline(args.min_align_score, color='k', linestyle='dashed')
    ax.set(xlabel='Alignment Score', ylabel='Observations')
    plt.savefig(filename)

plot_alignment_histogram(output['scores'].groupby(level='score').sum())

#merge_with_metadata_file(pd.concat([tallies, info], axis=1))

totals = tallies.sum()
reads = totals['Initial']
percents = defaultdict(float, totals/reads)

Log("""
Summary of the {:.2f}M processed reads of {:}:
-----------------------------------------------------------------------------
{Clustered:.1%} of reads will be used.
{Unaligned:.1%} did not align well to the master read and {PhiX:.1%} was PhiX.
{Wrong Barcode Length:.1%} had an inappropriate barcode length, while {Residual N:.1%} had an unfixable 'N' base. 
""".format(reads*1e-6, args.base_dir, **percents), print_line=True)

