#!/usr/bin/env python3

import argparse
from fastq import identify_RNA_barcode, csv_write
from ian_share import identify_variant_barcode

from shared_info import inerts, spike_BCs
import os
import pandas as p
p.options.mode.chained_assignment = None
import numpy as np
from collections import OrderedDict


outputs = [ 'full_clustered.csv', 'combined.csv']
abundance_metric = 'n10'

parser = argparse.ArgumentParser(description='Combine and superficially-analyze output of DADA2 Clustering.')
parser.add_argument('--full', dest='full', action='store_true')
parser.add_argument('-M', dest='multi_run', action='store_true', help='{:} already exists and has a 1st index column representing the multiple DADA2 runs that were combined.'.format(outputs[0]))
args = parser.parse_args()

def consolidate_clustered(directory='clustered/', ext='.csv', cluster_flank=7):
    Files = [f for f in os.listdir(directory) if ext in f]
    dfs = []
    keys = []
    matches = 0
    for f in Files:
        try:
            df = p.read_csv(directory+f, index_col=0)
        except: 
            print("Couldn't Load:", f)
            continue

        df.index = p.MultiIndex.from_tuples(df['sequence'].apply(identify_RNA_barcode, args=(cluster_flank,)))
        df.index.names = ['RNA', 'barcode']
        split = f.split('.csv')[0].split('_')
        if len(split) == 1:
            mouse, RNA = split[0], 'all'
        else:
            mouse, RNA = split if len(split) == 2 else split[::2]
        if RNA != 'all':
            matches += df.index.get_level_values('RNA').str.count(RNA).values.sum()
        dfs.append(df)
        keys.append((mouse, ))
    collected = p.concat(dfs, keys=keys, names=['Mouse'])
    if matches > 0:
        print("{:%} of sg IDs matched before and after clustering.".format(matches/len(collected)))
    csv_write(collected, outputs[0])
    return collected

MRB = ['Mouse', 'RNA', 'barcode']

def multi_run_similarity(df, threshold=0.95):
    print("The following Mice are found on multiple runs:")
    print("---------------------------------------")
    print('Mouse, Corr (Mean), Runs (reads) ...:')
    for mouse, df_m in df.groupby(level='Mouse'):
        S = df_m['n0'] + df_m['n1']
        M = S.unstack(level='run')
        M = M.fillna(0)
        if len(M.columns) == 1:
            continue
        C = M.corr()
        unraveled = C.values[np.triu_indices(len(C), 1)]
        print(mouse + ',', '{:.4f}'.format(unraveled.mean()), ',', ' '.join([run + ' ({:.1f}m)'.format(reads*1e-6) for run, reads in M.sum().items()]))

def prod(S): return np.multiply.reduce(S.values)

merge = dict(abundance=sum, n0=sum, n1=sum, nunq=sum, pval=prod, birth_pval=prod, birth_ham=min)

if args.multi_run:
    rMRB = ['run'] + MRB
    init = p.read_csv(outputs[0]+'.bz2', index_col=np.arange(4), dtype=dict(zip(rMRB, 4*[str])))
    gb = init.groupby(level=rMRB)
    df = gb.agg(merge)
    multi_run_similarity(df)
else:
    df = p.read_csv(outputs[0]+'.bz2', index_col=[0, 1, 2]) if os.path.isfile(outputs[0] + '.bz2') and not args.full else consolidate_clustered()

gb = df.groupby(level=MRB)
final = gb.agg(merge)
print("{:.4%} of barcodes were unique.".format(len(gb)/len(final)))

final['n10'] = final['n0'] + final.pop('n1')

#SAVE IAN HERE


abundance = final.loc[:, abundance_metric]
mice_groups = abundance.groupby(level='Mouse')

#### SPIKE Info
true_spike = abundance.loc[:, 'Spike', spike_BCs]
true_spike_percentage = true_spike/true_spike.groupby(level=['Mouse']).sum()*100
ranks = true_spike.groupby(level='Mouse').rank().unstack()
print(ranks.describe())

### Duplication info ###########
def duplicate_barcode(S):
    counts = S.groupby(level='barcode').count()
    nonzero = counts[counts > 1]
    duplicate_barcodes = frozenset(nonzero.index.get_level_values('barcode').values)
    return S.index.get_level_values('barcode').isin(duplicate_barcodes)

barcodes = final.index.get_level_values('barcode')
duplications = OrderedDict([('library', duplicate_barcode(abundance)),
                            ('mouse',   mice_groups.transform(duplicate_barcode).astype(bool)),
                            ('RNA',     abundance.groupby(level="RNA").transform(duplicate_barcode).astype(bool))])

final.loc[:, 'barcode duplicated in'] = ''
for name, S in duplications.items():
    final.loc[S, 'barcode duplicated in'] = name

final['max fraction'] = abundance.groupby(level=["RNA", 'barcode']).transform(lambda S: S/S.max())
csv_write(final, outputs[1])

#### Defunct

def delta_size(abundance):
    theory_abundance = abundance.loc[:, inerts].values
    if len(theory_abundance) <= 1: 
        abundance.loc[:] = np.nan
        return abundance
    gb = abundance.groupby(level='RNA')
    out = gb.apply(lambda S: S - np.percentile(theory_abundance, ((S.rank().values - 0.5)*(100/len(S)) ).tolist() ) )
    return out

#metrics = { 'Cells'         : spike_normalization,
#            '% of total'    : lambda S: S/S.sum(), 
#           'relative to Inert': lambda df: df/( df.loc[:, inerts].sum()/len(inerts) ),
#            'size delta'    : delta_size}

#final['relative to inert'] = mice_groups.apply(lambda df: df/df.loc[:, safe_inerts].mean())
#final['relative to inert sum'] = mice_groups.apply(lambda df: df/( df.loc[:, safe_inerts].sum()/len(safe_inerts) ))
#final['relative to cohort'] = abundance.groupby(level=['Mouse', 'RNA']).apply(lambda df: df/df.mean())

#for name, metric in metrics.items():
#    final[name] = mice_groups.transform(metric)

