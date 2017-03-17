import pandas as pd
import numpy as np
import os, argparse
from tuba_seq.fastq import fastqDF, NW_fit, repair_N, to_bytes
from tuba_seq.shared import logPrint
from tuba_seq.blast import sleuth_DNAs
import params
from rpy2.robjects.packages import importr
from rpy2.robjects.pandas2ri import py2ri, ri2py
import warnings
from memory_profiler import profile

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

###############################################################################

args = parser.parse_args()
os.chdir(args.base_dir)

dada2 = importr("dada2")
R_base = importr('base')

if args.parallel:
    from tuba_seq.pmap import pmap as map 

immediate_truncation = slice(max(0, len(params.head) - params.alignment_flank), 
                             min(len(params.master_read), len(params.head) + params.barcode_length + params.alignment_flank ))

ref = params.master_read[immediate_truncation].replace('.', 'N')

b_ref = to_bytes(ref)

training_flank = params.training_flank
cluster_flank = params.cluster_flank

def NW_fit_ref(seq):
    return NW_fit(to_bytes(seq), b_ref, training_flank)

_start, _stop, max_score = NW_fit(b_ref, b_ref, training_flank)
_start, _stop, min_score = NW_fit(to_bytes( ref.replace('A', 'T').replace('G', 'C') ), to_bytes( ref.replace('T', 'A').replace('C', 'G') ), training_flank)

X_score = np.linspace(min_score, max_score, num=HISTOGRAM_BINS)
Y_score = np.zeros_like(X_score[1:])

min_align_score = max_score*params.min_align_score_frac

os.makedirs(params.preprocessed_dir, exist_ok=True)
os.makedirs(params.training_dir, exist_ok=True)

input_dir = params.panda_seq['merge_dir'] if args.merged else params.original_dir

files = [f for f in os.listdir(input_dir) if params.fastq_handle in f]

Log = logPrint(verbose=args.verbose)
Log('Processing {:} samples found in {:}.'.format(len(files), os.path.join(args.base_dir, input_dir)), print_line=True)
summary = {}
unknown_DNAs = {}

#for filename in files:
@profile
def run(filename):
    basename = os.path.basename(filename)
    sample = basename.split('.')[0]
    start = fastqDF.from_file('original/'+filename).drop_Illumina_filter()
    info = start.fastq_info()
    
    df = start.co_slice(immediate_truncation).drop_abnormal_lengths()
    gb = df.groupby('DNA')
    DNAs = pd.Series(list(gb.groups.keys()))
    isPhiX = ri2py(dada2.isPhiX(py2ri(DNAs))).astype(bool)
    no_phiX = DNAs.loc[~isPhiX]
    NW_fits = pd.DataFrame(dict(zip(no_phiX, map(NW_fit_ref, no_phiX))), index=['start', 'stop', 'score']).T
    NW_fits['abundance'] = gb['DNA'].count()

    #Y_score += np.histogram(NW_fits['score'], weights=NW_fits['abundance'], bins=X_score)[0]
    NW_fits['aligned'] = NW_fits['score'] >= min_align_score
    reasonable_fits = NW_fits.query("aligned and abs(stop - start - {0.barcode_length}) <= {0.allowable_deviation}".format(params))
    reasonable_DNAs = frozenset(reasonable_fits.index.values.tolist())
    reasonable = gb.filter(lambda data: data.name in reasonable_DNAs )

    #heads
    #barcodes
    #tails
    #fake_header


    def training_trim(df):
        DNA = df['DNA'].values[0]
        start = reasonable_fits.loc[DNA, 'start']
        stop = reasonable_fits.loc[DNA, 'stop']
        left_slice = slice(start - training_flank, start)
        right_slice = slice(stop, stop + training_flank)
        df['DNA'] = DNA[left_slice]+DNA[right_slice]
        left_qc = df['QC'].str.slice(left_slice.start, left_slice.stop)
        right_qc = df['QC'].str.slice(right_slice.start, right_slice.stop)
        df['QC'] = left_qc.str.cat(right_qc) 
        return df

    def clustering_trim_and_repair(df):
        DNA = df['DNA'].values[0]
        start = reasonable_fits.loc[DNA, 'start']
        if 'N' in DNA:
            DNA = repair_N(to_bytes(DNA), b_ref)

        slicer = slice(start - cluster_flank, start + params.barcode_length + cluster_flank)
        df['DNA'] = DNA[slicer]
        df['QC'] = df['QC'].str.slice(slicer.start, slicer.stop)
        return df
   
    def trim_reads(gb, trim_function):          
        dfs = [df for _dna, df in gb]                                       # Using map function instead of apply for parallelization 
        trimed_df = fastqDF(pd.concat(list(map(trim_function, dfs))))       
        return trimed_df.drop_degenerate().drop_abnormal_lengths()

    def save_reads(df, directory):
        filename = os.path.join(directory, sample)
        df.write(filename, compression=None if args.derep else 'gzip')
        if args.derep:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore") 
                derep = dada2.derepFastq(filename+fastqDF.fastq_handle, verbose=args.verbose)
            R_base.saveRDS(derep, file=filename+'.RData')
            os.remove(filename+fastqDF.fastq_handle)

    next_gb = reasonable.groupby("DNA")
    
    save_reads(trim_reads(next_gb, training_trim), params.training_dir)
    cluster_df = trim_reads(next_gb, clustering_trim_and_repair) 
    EEs = cluster_df.expected_errors()
    cluster_df = cluster_df.select_reads(EEs <= params.maxEE)
    save_reads(cluster_df, params.preprocessed_dir) 

    init_reads = len(start)
    tally = {'Clustered Reads'  : len(cluster_df), 
            'Initial Reads'     : init_reads,
            'Aligned Reads'     : NW_fits[['aligned', 'abundance']].prod(axis=1).sum(),
        'Unfiltered PhiX Reads' : len(df) - NW_fits['abundance'].sum(),
            'Perfect Length'    : reasonable_fits.query('stop - start == {:}'.format(params.barcode_length))['abundance'].sum(),
            'Expected Errors'   : EEs.sum()}
    
    tally['Truncated Reads'] = tally['Aligned Reads'] - len(reasonable)
    tally = pd.Series(tally)
    
    percents = tally/init_reads
    Log('Sample {:}: {:.2}M Reads, {Aligned Reads:.1%} aligned, {Unfiltered PhiX Reads:.1%} PhiX, {Clustered Reads:.1%} Cluster-able.'.format(sample, init_reads*1e-6, **percents))
    
    info.pop('Flowcell ID'), info.pop('Flowcell Lane')

    summary[sample] = tally.append(info)
    
    NW_fits['fraction'] = NW_fits['abundance']/init_reads
    unknown_DNAs[sample] = NW_fits.query("not aligned and fraction >= @args.fraction")

for f in files:
    run(f)

unknown_DNAs = pd.concaat(unknown_DNAs, names=['Sample']) 
if args.search_blast and len(unknown_DNAs) > 0: 
    gb = unknown_DNAs.groupby(level='DNA')
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

summary = pd.DataFrame(summary).T.apply(pd.to_numeric, errors='ignore')

def plot_alignment_histogram(filename='alignment_histogram.pdf'):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    dx = X_score[1] - X_score[0]
    ax.bar(X_score[:-1], Y_score, width=dx)
    ax.axvline(min_align_score, color='k', linestyle='dashed')
    ax.set(xlabel='Alignment Score', ylabel='Observations')
    plt.savefig(filename)

plot_alignment_histogram()

if len(summary['Instrument'].value_counts()) > 1:
    Log("This run contains fastq files from two different Illumina machines. This is not recommended.", print_line=True)

def merge_with_metadata_file(summary):
    pass

merge_with_metadata_file(summary)

totals = summary.select_dtypes(exclude=[object]).sum()
reads = totals['Initial Reads']
totals['Unaligned Reads'] = reads - totals['Aligned Reads']

EE_rate = totals['Expected Errors']/(totals['Clustered Reads']*(params.barcode_length + 2*params.cluster_flank))
percents = totals/reads

Log("""
Summary for {:}:
-----------------------------------------------------------------------------
Summary of the {:.2}M processed reads:
{Clustered Reads:.1%} will be used.
{Unaligned Reads:.1%} of reads did not align well to the master read ({Unfiltered PhiX Reads:.1%} was PhiX, {Truncated Reads:.1%} were truncated.)
The remainder generally have unfixable 'N' nucleotides or poor Phred Scores.  
The anticipated Phred-based Error Rate is {EE_rate:.4%} 
(Error training will undoubtedly find a larger error rate.)""".format(  args.base_dir, 
                                                        reads*1e-6, 
                                                        EE_rate=EE_rate,
                                                        **percents.to_dict()), print_line=True)
Log.close()

