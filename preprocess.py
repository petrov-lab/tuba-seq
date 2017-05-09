import pandas as pd
import os, argparse
from tuba_seq.fastq import fastqDF, repair_N, find_start_stop, nw
from tuba_seq.shared import logPrint
from tuba_seq.blast import sleuth_DNAs
import params
from rpy2.robjects.packages import importr
from collections import defaultdict
from Bio.Seq import Seq
import warnings

############### Input Parameters ########################

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.")
parser.add_argument('-b', "--base_dir", type=str, help='Base directory to work within. This directory must contain a folder entitled "{:}" containing all FASTQ files to process.'.format(params.original_dir), default=os.getcwd())
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multithreaded operation')
parser.add_argument('-s', '--search_blast', action='store_true', help='Use NCBI BLAST algorithm to identify contaminations in samples')
parser.add_argument('-l', '--local_blast', action='store_true', help='Use local NCBI BLAST+ algorithm to accelerate searching (if present), see tuba_seq/blast.py.')
parser.add_argument('-d', '--derep', action='store_true', help='De-replicate output fastQ files for DADA2 to minimize file sizes.')
parser.add_argument('-f', '--fraction', type=float, default=0.05, help='Minimum fraction of total reads to elicit a BLAST-search of an unknown sequence.')
parser.add_argument('-m', '--merged', action='store_true', help='Use merged reads in panda_seq.merge_dir.')
HISTOGRAM_BINS = 30
COMPRESSION = 'bz2'
###############################################################################

args = parser.parse_args()
os.chdir(args.base_dir)

dada2 = importr("dada2")
R_base = importr('base')

fastq_ext = params.fastq_handle
if not args.derep:
    if COMPRESSION == 'gzip':
        from gzip import open       # We don't need to specify binary mode...it will automatically assume.    
        fastq_ext += '.gz'
    elif COMPRESSION == 'bz2':
        from bz2 import open
        fastq_ext += '.bz2'

if args.parallel:
    from tuba_seq.pmap import low_memory_pmap as map 

ref = params.master_read.replace('.', 'N')

c_ref = ref.encode('ascii')

training_flank = params.training_flank
cluster_flank = params.cluster_flank
alignment_flank = params.alignment_flank

c_ref_scoring = c_ref[len(params.head) - alignment_flank:len(params.head) + params.barcode_length + alignment_flank]

max_score = nw.char_score(c_ref_scoring, c_ref_scoring)
opposite_ref_scoring = str(Seq(c_ref_scoring.decode('ascii')).complement())
min_score = nw.char_score(c_ref_scoring, opposite_ref_scoring.encode('ascii'))

min_align_score = max_score*params.min_align_score_frac

os.makedirs(params.preprocessed_dir, exist_ok=True)
os.makedirs(params.training_dir, exist_ok=True)

input_dir = params.merge_params['merge_dir'] if args.merged else params.original_dir

files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if params.fastq_handle in f]

def process_read_set(QCs, read_constructor, training_file, cluster_file):
    dna = QCs.index.values[0]
    c_dna = dna.encode('ascii')
    start, stop = find_start_stop(c_dna, c_ref)
    c_dna_scoring = c_dna[start - params.alignment_flank:stop + params.alignment_flank]
    score = nw.char_score(c_dna_scoring, c_ref_scoring)
    if args.merged and params.compress_PHRED:
        QCs = QCs.str.translate(params.PHRED_compressor)
    
    if score < min_align_score:
        return 'PhiX' if dada2.isPhiX(c_dna)[0] else 'Unaligned', score
    if abs((stop - start) - params.barcode_length) > params.allowable_deviation:
        return 'Wrong Barcode Length', score
    training_DNA = dna[start - training_flank:start]+dna[stop:stop + training_flank]
    if 'N' not in training_DNA and len(training_DNA) == 2*training_flank:
        training_QCs =  QCs.str.slice(start - training_flank, start).str.cat(
                        QCs.str.slice(stop, stop + training_flank))
        training_file.write(read_constructor(training_DNA, training_QCs)) 
    if 'N' in dna:
        dna = repair_N(c_dna, c_ref)
    cluster_DNA = dna[start - cluster_flank:start + params.barcode_length + cluster_flank]
    if 'N' not in cluster_DNA and len(cluster_DNA) == (params.barcode_length + 2*cluster_flank):
        cluster_file.write(read_constructor(cluster_DNA, 
                                            QCs.str.slice(start - cluster_flank, start + params.barcode_length + cluster_flank)))
        return 'Clustered', score
    else:
        return 'Residual N' if 'N' in cluster_DNA else 'Insufficient Flank', score

Log = logPrint(verbose=args.verbose)
Log('Processing {:} samples found in {:}.'.format(len(files), os.path.join(args.base_dir, input_dir)), print_line=True)

def process_sample(filename):
    reads = fastqDF.from_file(filename, fake_header=True, use_Illumina_filter=True)
    info = reads.info
    sample = info['Sample']
    init_reads = len(reads)
    gb = reads.set_index("DNA")['QC'].groupby(level='DNA') 
    filenames = [os.path.join(dir, sample)+fastq_ext for dir in (params.training_dir, params.preprocessed_dir)]
    with open(filenames[0], 'w') as training_file, open(filenames[1], 'w') as cluster_file:
        reads_tuples = gb.agg(process_read_set, reads.construct_read_set, training_file, cluster_file)
    reads_ix = pd.MultiIndex.from_tuples(reads_tuples, names=['outcome', 'score'])
    counts = gb.count().reset_index()
    counts.index = reads_ix
    counts = counts.set_index('DNA', append=True)['QC']
    counts.name = 'Reads' 
    tallies = counts.groupby(level='outcome').sum()
    percents = defaultdict(float, tallies/init_reads)
    Log('Sample {:}: {:.2f}M Reads, {Unaligned:.1%} unaligned, {PhiX:.1%} PhiX, {Wrong Barcode Length:.1%} inappropriate barcode length, {Clustered:.1%} Cluster-able.'.format(sample, init_reads*1e-6, **percents))
    tallies['Initial'] = init_reads
    return dict(tallies=tallies,
                info=info, 
                unknown_DNAs=counts.reset_index('outcome').query("outcome == 'Unaligned' and Reads > {:}".format(int(args.fraction*init_reads)))['Reads']/init_reads,
                scores=counts.groupby(level='score').sum())

_reports = map(process_sample, files)

dfs = {r['info'].loc['Sample']:r for r in _reports}
meta_datas = ['tallies', 'unknown_DNAs', 'info', 'scores']

output = {md:pd.concat({d['info'].loc['Sample']:d[md] for d in _reports}, names=['Sample']) for md in meta_datas}

def derep(sample, directories=(params.training_dir, params.preprocessed_dir)):
    filenames = [os.path.join(Dir, sample+fastq_ext) for Dir in directories]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") 
        for f in filenames:
            try:
                derep = dada2.derepFastq(f, verbose=args.verbose)
                R_base.saveRDS(derep, file=f.replace(fastq_ext, '.rds'))
            except Exception: 
                Log("Could not derep "+f)
            else:
                os.remove(f)

if args.derep:
    print("De-replicating outputs for DADA2...")
    for _ in map(derep, dfs.keys()):
        pass

if args.search_blast and len(output['unknown_DNAs']) > 0: 
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

if len(info['Instrument'].value_counts()) > 1:
    Log("This run contains fastq files from two different Illumina machines. This is not recommended.", print_line=True)

def merge_with_metadata_file(summary):
    pass

def plot_alignment_histogram(data, filename='alignment_histogram.pdf', hist_kwargs={'bins':30}):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    ax.hist(data.index.values, weights=data.values, **hist_kwargs)
    ax.axvline(min_align_score, color='k', linestyle='dashed')
    ax.set(xlabel='Alignment Score', ylabel='Observations')
    plt.savefig(filename)

plot_alignment_histogram(output['scores'].groupby(level='score').sum())

merge_with_metadata_file(pd.concat([tallies, info], axis=1))

totals = tallies.sum()
reads = totals['Initial']
percents = defaultdict(float, totals/reads)

Log("""
Summary of the {:.2f}M processed reads of {:}:
-----------------------------------------------------------------------------
{Clustered:.1%} of reads will be used.
{Unaligned:.1%} did not align well to the master read ({PhiX:.1%} was PhiX.)
{Wrong Barcode Length:.1%} had an inappropriate barcode length, while {Residual N:.1%} had an unfixable 'N' base. 
""".format(reads*1e-6, args.base_dir, **percents), print_line=True)
Log.close()

