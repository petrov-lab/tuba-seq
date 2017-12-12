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
from matplotlib import pyplot as plt

fastq_ext = '.fastq'
histogram_filename = 'alignment_histogram.pdf'

default_master_read = ''.join((
    'GCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGA',# forward flank
    '........',                                         # sgID
    'AA.....TT.....AA.....',                            # random barcode
    'ATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAAGAACCT')) # aft flank

############### Input Parameters ########################

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('input_dir', type=str, help='Input directory with fastq files.') 
parser.add_argument('--master_read', type=str, default=default_master_read, 
    help="Outline of the amplicon sequence, degenerate bases can be either 'N' or '.'")

parser.add_argument('-t', '--training_dir', default='training/', help='Directory to save files for error training.') 
parser.add_argument('-o', '--output_dir', default='preprocessed/', help='Directory to save files for barcode clustering.') 
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
parser.add_argument('--training_flank', type=int, default=18, help='# of bases flanking the degenerate region to be used to develop the DADA2 error model.')
parser.add_argument('-M', '--min_align_score', type=float, default=0.6, help='Minimum alignment score needed to keep read, Range [0, 1).')
parser.add_argument('--compression', default='bz2', choices=['bz2', 'gzip', 'none'], help='Compression algorithm for saved file.')
###############################################################################
args = parser.parse_args()

master_read = MasterRead(args.master_read, args)
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
else:
    raise ValueError("Unknown compression: "+args.compression)

if args.parallel:
    from tuba_seq.pmap import low_memory_pmap as map

for Dir in output_dirs:
    os.makedirs(Dir, exist_ok=True)

files = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if fastq_ext in f]
max_file_size = 1e9 if args.parallel else 2e9

Log = logPrint(args)
Log('Processing {:} samples found in {:}.'.format(len(files), args.input_dir), print_line=True)

manager = Manager()
master_read.scores = manager.dict()
master_read.unaligned = manager.dict()
master_read.alignments = manager.dict()
master_read.instruments = manager.dict()


def process_fastq(filename):
    sample = os.path.basename(filename.partition(fastq_ext)[0])
    filenames = [os.path.join(Dir, sample)+fastq_ext+compression_ext for Dir in output_dirs]
    if args.skip and all([os.path.isfile(f) for f in filenames]):
        Log("Skipping "+sample+" (already exists).")
        return 

    if filename[-3:] == '.gz' or filename[-5:] == '.gzip':
        from gzip import open as read_open
    elif filename[-4:] == '.bz2':
        from bz2 import open as read_open
    elif filename.partition(fastq_ext)[2] == '': 
        from builtins import open as read_open
    else:
        raise ValueError(filename + " is not a valid file extension to open.")
    
    with read_open(filename) as input_file:
        outcomes = master_read.iter_fastq(sample, filenames, input_file)
    reads = outcomes.sum()
    Log('Sample {:} ({:.2f}M Reads): '.format(sample, reads*1e-6)+
        ','.join(['{:.1%} {:}'.format(num/reads, name) for name, num in outcomes.iteritems() if num > 0])+'.')

    if args.derep:
        Log("De-replicating outputs for "+sample+" ...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") 
            for f in filenames:
                try:
                    derep = dada2.derepFastq(f, verbose=args.verbose)
                    R_base.saveRDS(derep, file=f.replace(fastq_ext, '.rds'))
                except Exception: 
                    print("Could not derep "+f)
                else:
                    os.remove(f)
    
    return sample, outcomes

outcomes = pd.DataFrame(dict(map(process_fastq, files))).T

scores = pd.DataFrame(dict(master_read.scores)).sum(axis=1) 
unaligned = pd.Series(dict(master_read.unaligned))

total_reads = outcomes.sum().sum()

PhiX = pandas2ri.ri2py(dada2.isPhiX(pandas2ri.py2ri(unaligned.index))) == 1
unknown_DNAs = (unaligned.loc[~PhiX]/total_reads).loc[lambda x: x >= args.fraction]

if args.search_blast and len(unknown_DNAs) > 0: 
    from tuba_seq.blast import sleuth_DNAs
    DNAs = unknown_DNAs.index.values
    Log("BLAST-searching", len(DNAs), "common unknown sequence(s)...", True)
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

ax = plt.gca()
ax.hist(scores.index.values, weights=scores.values, bins=numpy.linspace(0, 1, 40))
ax.axvline(args.min_align_score, color='k', linestyle='dashed')
ax.set(xlabel='Alignment Score', ylabel='Observations')
plt.savefig(histogram_filename)

totals = outcomes.sum()
totals['PhiX'] = unaligned.loc[PhiX].sum()

Log("Summary of the {:.2f}M processed reads in {:}:".format(total_reads*1e-6, args.input_dir), True, header=True)
Log((totals/total_reads).to_string(float_format='{:.2%}'.format), True)

