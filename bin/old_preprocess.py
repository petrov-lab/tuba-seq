#!/usr/bin/env python3
import pandas as pd
import os, numpy, argparse
from multiprocessing import Manager
from tuba_seq.fastq import MasterRead
from tuba_seq.shared import logPrint
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import warnings

fastq_ext = '.fastq'
histogram_filename = 'alignment_histogram.pdf'

default_master_read = ''.join((
    'GCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA',# forward flank
    '........',                                         # sgID
    'AA.....TT.....AA.....',                            # random barcode
    'ATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCT')) # aft flank

############################ Input Parameters #################################
parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('input_dir', type=str, help='Input directory with fastq files.') 
parser.add_argument('--single', action='store_true', help='Analyze single FASTQ file; The positional argument should be a *file*, not a directory.')
parser.add_argument('--master_read', type=str, default=default_master_read, 
    help="Outline of the amplicon sequence, degenerate bases can be either 'N' or '.'. --trim and --symmetric_flanks depend on this being the full length FASTQ sequence that you expect after merging reads.")
parser.add_argument('-t', '--training_dir', default='training', help='Directory to save files for error training.') 
parser.add_argument('-o', '--output_dir', default='preprocessed', help='Directory to save files for barcode clustering.') 
parser.add_argument('--tally_filename', default='mutation_tallies.csv', help='CSV file of mutation tallies.')

#### NEED TO BE *PER SAMPLE*!
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multi-threaded operation')
parser.add_argument('-s', '--search_blast', action='store_true', help='Use NCBI BLAST algorithm to identify contaminations in samples')
parser.add_argument('-l', '--local_blast', action='store_true', 
    help='Use local NCBI BLAST+ algorithm to accelerate searching (if present), see tuba_seq/blast.py.')
parser.add_argument('-d', '--derep', action='store_true', help='De-replicate output fastQ files for DADA2 to minimize file sizes.')
parser.add_argument('-f', '--fraction', type=float, default=0.01, help='Minimum fraction of total reads to elicit a BLAST-search of an unknown sequence.')
parser.add_argument('-k',  '--skip', action='store_true', help='Skip files that already exist in output directories.')
parser.add_argument('-a', '--allowable_deviation', type=int, default=4, help="Length of Indel to tolerate before discarding reads.")
parser.add_argument('--alignment_flank', type=int, default=22, help='# of bases flanking the degenerate region to be used to score the quality of the read.')
parser.add_argument('--training_flank', type=int, default=22, help='# of bases flanking the degenerate region to be used to develop the DADA2 error model.')
parser.add_argument('-M', '--min_align_score', type=float, default=0.6, help='Minimum alignment score needed to keep read, Range [0, 1).')
parser.add_argument('--compression', default='bz2', choices=['bz2', 'gz', 'lzma', 'none'], help='Compression algorithm for saved file.')
parser.add_argument('--symmetric_flanks', action='store_true', help='Automatically symmetrizes the length of the forward and aft flanks to the barcode--saves memory/time.')
parser.add_argument('--trim', type=int, default=5, help='Trim reads by specified length on both flanks--works with --symmetric_flanks.')
parser.add_argument('--max_GB', type=int, default=32, help='Maximum quantity of memory to use for alignment memoization (in Gigabytes).') 
###############################################################################
args = parser.parse_args()

master_read = MasterRead(args.master_read, args)
output_dirs = args.training_dir, args.output_dir
dada2 = importr("dada2")
R_base = importr('base')

compression = '' if args.derep or args.compression == 'none' else '.'+args.compression

if args.parallel:
    from tuba_seq.pmap import pmap as map

input_fastqs = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if fastq_ext in f] if not args.single else [args.input_dir]

Log = logPrint(args)
if not args.single:
    Log('Processing {:} samples found in {:}.'.format(len(input_fastqs), args.input_dir), print_line=True)

manager = Manager()
master_read.scores = manager.dict()
if args.search_blast:
    master_read.unaligned = manager.dict()
master_read.alignments = manager.dict()

BYTES_PER_ENTRY = 100
GIGABYTE = pow(2, 30)
master_read.max_alignments = args.max_GB/BYTES_PER_ENTRY*GIGABYTE
master_read.instruments = manager.dict()
master_read.PHRED_tally = manager.dict()

def process_fastq(input_fastq):
    sample = os.path.basename(input_fastq.partition(fastq_ext)[0])
    filenames = [os.path.join(Dir, sample)+fastq_ext+compression for Dir in output_dirs]
    if args.skip and all([os.path.isfile(f) for f in filenames]):
        Log("Skipping "+sample+" (already exists).")
        return sample, pd.Series(dict(Filtered=0, Unaligned=0, Wrong_Barcode_Length=0, Residual_N=0, Insufficient_Flank=0, Clustered=0))
    
    outcomes = master_read.iter_fastq(sample, filenames, input_fastq)
    reads = outcomes.sum()
    Log('Sample {:} ({:.2f}M Reads): '.format(sample, reads*1e-6)+
        ','.join(['{:.1%} {:}'.format(num/reads, name) for name, num in outcomes.iteritems() if num > 0])+'.')
    if outcomes['Clustered'] == 0:
        Log('There were no passable reads in {:}. Deleting output files...', True)
        list(map(os.remove, filenames))
    return sample, outcomes

outcomes = pd.DataFrame(dict(map(process_fastq, input_fastqs))).T

scores = pd.DataFrame(dict(master_read.scores)).sum(axis=1) 
total_reads = outcomes.sum().sum()

if len(master_read.alignments) == master_read.max_alignments:
    Log('Memoization of read alignments maxed-out; increase `max_GB` for better performance.', True)

PHRED_tally = pd.Series(dict(master_read.PHRED_tally))
PHRED_tally.index.names = ['mutation', 'PHRED']
PHRED_tally.unstack('PHRED').fillna(0).astype(int).to_csv(args.tally_filename)

if args.search_blast:
    unaligned = pd.Series(dict(master_read.unaligned))
    unaligned.index = unaligned.index.astype(str)
    PhiX = pandas2ri.ri2py(dada2.isPhiX(pandas2ri.py2ri(unaligned.index))) == 1
    non_PhiX = unaligned.loc[~PhiX]
    unknown_DNAs = (non_PhiX/total_reads).loc[lambda x: x >= args.fraction]

    if len(unknown_DNAs) > 0: 
        from tuba_seq.blast import sleuth_DNAs
        DNAs = unknown_DNAs.index.values
        Log("BLAST-searching {:} common unknown sequence(s)...".format(len(DNAs)), True)
        searches = sleuth_DNAs(DNAs, local_blast=args.local_blast)
        if len(searches) == 0:
            Log("Could not find any acceptable matches.", True)
        else:
            Log("                -- Best Match --                          | E-score | % of Reads", True)
            for dna, row in searches.iterrows():
                Log("{short_title:<60} {e:1.0e}    {:.2%}".format(unknown_DNAs[dna], **row),True)

instruments = pd.Series(dict(master_read.instruments))
if len(instruments.value_counts()) > 1:
    Log("This run contains fastq files from two different Illumina machines. This is not recommended.", True)
    Log(instruments, True)

try:
    from tuba_seq.graphs import plt
    ax = plt.gca()
    ax.hist(scores.index.values, weights=scores.values, bins=numpy.linspace(0, 1, 40))
    ax.axvline(args.min_align_score, color='k', linestyle='dashed')
    ax.set(xlabel='Alignment Score', ylabel='Observations')
    plt.savefig(histogram_filename)
except Exception as e:
    print("Couldn't create alignment score histogram, perhaps you need to configure the matplotlib backend?")
    print(e)

totals = outcomes.sum()
if args.search_blast:
    totals['PhiX (subset of Unaligned)'] = unaligned.loc[PhiX].sum()

Log("Summary of the {:.2f}M processed reads in {:}:".format(total_reads*1e-6, args.input_dir), True, header=True)
Log((totals/total_reads).to_string(float_format='{:.2%}'.format), True)

if args.derep:
    Log("De-replicating outputs...")
    filenames = [os.path.join(Dir, sample)+fastq_ext+compression for Dir in output_dirs for sample in outcomes.index.values]
    sizes = list(map(os.path.getsize, filenames))
    files_by_size = pd.Series(sizes, index=filenames).sort_values() 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        for f in files_by_size.keys():
            try:
                print(f)
                derep = dada2.derepFastq(f, verbose=args.verbose)
                R_base.saveRDS(derep, file=f.replace(fastq_ext, '.rds'))
            except Exception as e: 
                Log("Could not derep {:}:\n{:}".format(f, e), True)
            #else:
                #os.remove(f)
    
