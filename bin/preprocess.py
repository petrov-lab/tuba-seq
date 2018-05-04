#!/usr/bin/env python3
import pandas as pd
import os, numpy, argparse, sys, warnings
from tuba_seq.fastq import MasterRead
from tuba_seq.shared import logPrint, smart_open
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

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
parser.add_argument('--master_read', type=str, default=default_master_read, 
    help="Outline of the amplicon sequence, degenerate bases can be either 'N' or '.'. --trim and --symmetric_flanks depend on this being the full length FASTQ sequence that you expect after merging reads.")
parser.add_argument('-t', '--training_dir', default='training', help='Directory to save files for error training.') 
parser.add_argument('-o', '--output_dir', default='preprocessed', help='Directory to save files for barcode clustering.') 
parser.add_argument('-u', '--unaligned_dir', default='unaligned', help='Directory to save unaligned reads')
parser.add_argument("-v", "--verbose", help='Output more', action="store_true")
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
parser.add_argument('-M', '--min_align_score', type=float, default=0.65, help='Minimum alignment score needed to keep read, Range [0, 1).')
parser.add_argument('--compression', default='bz2', choices=['bz2', 'gz', 'lzma', 'none'], help='Compression algorithm for saved file.')
parser.add_argument('--trim', default='symmetric', help='Nucleotides to immediately trim from the amplicon reads before searching for the barcode--trimming accelerates runtime. Can be two integers--a start and stop position, `none`, or `symmetric`, which truncates the read such that the barcode is exactly in the middle of the read.')
###############################################################################
args = parser.parse_args()
Log = logPrint(args)

master_read = MasterRead(args.master_read, args)

dada2 = importr("dada2")
R_base = importr('base')

single_file = '.fastq' in args.input_dir

compression = '' if args.derep or args.compression == 'none' else '.'+args.compression

if args.parallel and not single_file:
    from tuba_seq.pmap import pmap as map

input_fastqs = [os.path.join(args.input_dir, f) for f in os.listdir(args.input_dir) if fastq_ext in f] if not single_file else [args.input_dir]

def get_instrument(filename):
    with smart_open(filename) as f: 
        return f.readline().decode('ascii').split(':')[0]

instruments = pd.Series({input_fastq:get_instrument(input_fastq) for input_fastq in input_fastqs})

if len(instruments.value_counts()) > 1:
    Log("{input_dir} contains fastq files from different Illumina machines. This is not recommended.", True)
    Log(instruments, True)
    sys.exit()

if not single_file:
    Log('Processing {:} samples found in {:}.'.format(len(input_fastqs), args.input_dir), True)

def process_fastq(ix):
    sample = samples[ix]
    input_fastq = input_fastqs[ix]
    fastqs = fastq_outputs[ix]
    output_files = fastqs + fasta_outputs[ix:ix+1]
    if args.skip and all([os.path.isfile(f) for f in fastqs]):
        Log(sample+" already exists, skipping...")
        return 
    
    if args.parallel and single_file:
        from tuba_seq.pmap import fastq_map_sum
        outcomes, scores, bad_lengths = fastq_map_sum(input_fastq, output_files, master_read.iter_fastq)
    else:
        from tuba_seq.fastq import IterFASTQ
        outcomes, scores, bad_lengths = master_read.iter_fastq(IterFASTQ(input_fastq), output_files)
    reads = outcomes.sum()
    Log('Sample {:} ({:.2f}M Reads): '.format(sample, reads*1e-6)+
        ','.join(['{:.1%} {:}'.format(num/reads, name) for name, num in outcomes.iteritems() if num > 0])+'.')
    if outcomes['Clustered'] == 0:
        Log('There were no passable reads in {:}. Deleting output files...'.format(input_fastq), True)
        list(map(os.remove, output_files))
        return 

    if args.derep:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") 
            for f in fastqs:
                try:
                    print(f)
                    derep = dada2.derepFastq(f, verbose=args.verbose)
                    R_base.saveRDS(derep, file=f.replace(fastq_ext, '.rds'))
                except Exception as e: 
                    Log("Could not derep {:}:\n{:}".format(f, e), True)
                else:
                 os.remove(f)
    return outcomes, scores, bad_lengths

samples = [os.path.basename(input_fastq.partition(fastq_ext)[0]) for input_fastq in input_fastqs]
fastq_outputs = [[os.path.join(Dir, sample+fastq_ext+compression) for Dir in [args.training_dir, args.output_dir]] for sample in samples]
fasta_outputs = [os.path.join(args.unaligned_dir, sample+'.fasta'+compression) for sample in samples]
outputs = list(filter(lambda output: output is not None, map(process_fastq, range(len(samples)))))

if not outputs:
    Log("No files were processed.")
    sys.exit()

outcome_totals, score_totals, bad_barcode_length_totals = [sum(output_set) for output_set in zip(*outputs)]
total_reads = outcome_totals.sum()

if args.search_blast:  
    from collections import Counter
    unaligned = pd.Series(sum([Counter(open(fasta, 'rb')) for fasta in fasta_outputs], Counter()))
    unaligned.index = unaligned.index.astype(str).str.slice(0, -1)

    PhiX = pandas2ri.ri2py(dada2.isPhiX(pandas2ri.py2ri(unaligned.index))) == 1
    outcome_totals['PhiX (subset of Unaligned)'] = unaligned.loc[PhiX].sum()
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

try:
    from tuba_seq.graphs import plt, text_color_legend
    ax = plt.gca()
    bar_width = score_totals.index.values[1]
    removed = score_totals[score_totals.index < args.min_align_score]
    passed = score_totals[score_totals.index >= args.min_align_score]
    ax.bar(passed.index.values, passed.values, bar_width, label='Passed Alignment Filter')
    ax.bar(removed.index.values, removed.values, bar_width, label='Failed Alignment Filter')
    text_color_legend(ax, bbox_to_anchor=(0, 1), loc='upper left')
    ax.set(xlabel='Alignment Score', ylabel=score_totals.name, xlim=[0, 1+bar_width/2])
    plt.savefig(histogram_filename)
except Exception as e:
    print("Couldn't create alignment score histogram, perhaps you need to configure the matplotlib backend?")
    print(e)


Log("Summary of the {:.2f}M processed reads in {:}:".format(total_reads*1e-6, args.input_dir), True, header=True)
Log((outcome_totals/total_reads).to_string(float_format='{:.2%}'.format), True)

bad_lengths = bad_barcode_length_totals.sum()
if bad_lengths > 0:
    Log("Most common bad barcode length: {:} ({:.0%} of all bad lengths)".format(bad_barcode_length_totals.idxmax(), bad_barcode_length_totals.max()/bad_lengths), True)



