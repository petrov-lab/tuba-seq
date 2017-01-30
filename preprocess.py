import pandas as pd
import numpy as np
import os, argparse
from fastq import fastqDF, NW_fit, repair_N
from datetime import datetime
from matplotlib import pyplot as plt
from array import array
import params

USE_PMAP = True

start = datetime.now()

master_read = params.master_read
training_flank = params.training_flank

head = master_read[:master_read.find('.')]
tail = master_read[master_read.rfind('.')+1:]
degen_len = len(master_read) - len(head) - len(tail)

head_tail_diff = len(head) - len(tail)                      #  27 
symm_head_tail_trunc = slice(head_tail_diff, len(master_read))

ref = master_read[symm_head_tail_trunc].replace('.', 'N')

def to_bytes(s): 
    return array('B', s.encode('ascii'))

b_ref = to_bytes(ref)

_start, _stop, max_score = NW_fit(b_ref, b_ref, training_flank)
_start, _stop, min_score = NW_fit(to_bytes( ref.replace('A', 'T').replace('G', 'C') ), to_bytes( ref.replace('T', 'A').replace('C', 'G') ), training_flank)
X_score = np.linspace(min_score, max_score, num=params.hist_number)
Y_score = np.zeros_like(X_score[1:])

min_align_score = max_score*params.min_align_score_frac

parser = argparse.ArgumentParser(description="Prepare FASTQ files for DADA training & clustering.")
parser.add_argument("--fastq_dir", type=str, help='Directory containing FASTQ files in "original/" folder.', default=os.getcwd())
args = parser.parse_args()
root_dir = args.fastq_dir.rstrip('/')
os.chdir(root_dir)

os.makedirs(params.preprocessed, exist_ok=True)
os.makedirs(params.training, exist_ok=True)

files = [f for f in os.listdir('original/') if params.fastq_handle in f]

print('Processing {:} Files in {:}/original/.'.format(len(files), root_dir))

def NW_fit_ref(seq):
    return NW_fit(to_bytes(seq), b_ref, training_flank)

if USE_PMAP:
    from multiprocessing import Pool, cpu_count
    threads = cpu_count() - 1
    map = Pool(processes=threads).map

summary = {}

for f in files:
    short_f = f.split(params.fastq_handle)[0]
    df = fastqDF('original/'+f, symm_head_tail_trunc)
    df = df.join(pd.DataFrame(map(NW_fit_ref, df['DNA'].values), columns=['degen start', 'degen stop', 'score'], dtype=np.int16))
    df.__class__ = fastqDF 
    
    Y_score += np.histogram(df['score'].values, bins=X_score)[0]
    
    df['good_align'] = df['score'] > min_align_score 
    df['proper_degen'] = (df['degen stop'] - df['degen start'] - degen_len).abs() <= params.allowable_deviation
    
    reasonable = df.loc[df[['good_align', 'proper_degen']].all(axis=1), :]
    reasonable.__class__ = fastqDF
    counts = pd.Series(dict(reads=len(df), aligned=df['good_align'].sum(), median_score=df.loc[df['good_align'], 'score'].median()))
    counts['truncated'] = counts['aligned'] - len(reasonable)

    train_starts = reasonable['degen start'].apply(lambda start: slice(start-training_flank, start))
    train_stops  = reasonable['degen stop'].apply(lambda stop:  slice(stop, stop + training_flank))
    train = reasonable.vslice(train_starts, train_stops, enforce=True).dropDegen()
    train.write(params.training+f)

    need_repair_ix = reasonable['DNA'].str.contains("N")
    need_repair = reasonable.loc[need_repair_ix]
    reasonable.loc[need_repair_ix, 'DNA'] = need_repair['DNA'].apply(to_bytes).apply(repair_N, args=(b_ref, )) 

    cluster_slices = reasonable['degen start'].apply(lambda start: slice(start - params.cluster_flank, start + degen_len + params.cluster_flank))
    sliced = reasonable.vslice(cluster_slices, enforce=True)
    cluster = sliced.dropDegen()
    counts['unrepaired'] = len(sliced) - len(cluster)

    clean, EEs = cluster.expected_errors(maxEE=params.maxEE)
    counts['EEs'] = int(round(EEs))
    
    counts['High_EE'] = len(cluster) - len(clean) 

    final = clean.drop_abnormal_lengths()
    counts['final'] = len(final)
    final.write(params.preprocessed+('{short_f}.fastq'.format(**locals())))
    summary[short_f] = counts

summary = pd.DataFrame(summary).T
summary.to_csv('preprocessing_stats.csv')

runtime = datetime.now() - start
totals = summary.sum()

reads = totals['reads']
EE_rate = totals['EEs']/(totals['final']+totals['High_EE'])/(degen_len + 2*params.cluster_flank)

percents = totals/reads

ax = plt.gca()
ax.bar(X_score[:-1], Y_score, width=X_score[1] - X_score[0]) 
ax.axvline(min_align_score, color='k')
ax.set(xlabel='Alignment Score', ylabel='Observations')
plt.savefig('alignment_histogram.pdf') #, format='pdf', transparent=True, bbox_inches='tight')

print("""
Summary for {:} (took {:} h:m:s):
-----------------------------------------------------------------------------
Processed {:.2}M Reads.
{aligned:.1%} aligned to the Expected Sequence.
The median score was {relative_score:.1%} of Max.
{truncated:.1%} had an unreasonable-lengthed barcode.
{unrepaired:.1%} had unfixable N nucleotides.
{final:.1%} made it through. 
Their average Error Rate was {EE_rate:.4%}.""".format(  root_dir, 
                                                        str(runtime).split('.')[0], 
                                                        reads*1e-6, 
                                                        EE_rate=EE_rate,
                                                        relative_score=summary['median_score'].mean()/max_score,
                                                        **percents.to_dict()))
