#!/usr/bin/env python3
import argparse, os
import pandas as pd
import numpy as np
from tuba_seq.shared import logPrint
from pathlib import Path

tuba_seq_dir = os.path.dirname(__file__)

############################## Input #########################################

parser = argparse.ArgumentParser(   description="Consolidates mutation tallies directly measured from preprocess.py  & reports statistics.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('input_dir', type=Path, help='Input directory with mutation tally files.') 
parser.add_argument('--outfile', type=str, default='fullTransitionMatrix.csv', help='Directory to output full transition matrix.')

args = parser.parse_args()
Log = logPrint(args)

correct = ['A2A', 'C2C', 'G2G', 'T2T']
nucleotides = list('ACGT')

def reindex_and_clean(df):
    df = df.loc[(nucleotides, nucleotides), :]
    df.index = df.index.map('{0[0]}2{0[1]}'.format)
    df.index.names = ['Mutation']
    return df

data = pd.concat({f.stem:reindex_and_clean(pd.read_csv(f, index_col=[0, 1])) for f in args.input_dir.iterdir() if f.suffixes[0] == '.csv'}, names=['Sample'])
data.columns = data.columns.astype(int)
relevant = data.loc[:, data.sum() > 0]

full_trans = relevant.groupby(level='Mutation').sum()

# Now the report
gb = relevant.groupby(level='Sample')

def error_rate(col):
    return 1 - col.loc[:, correct].sum()/col.sum()

stats = gb.agg([error_rate, sum]).stack(level=0).sort_index()
stats.index.names = ['Sample', 'PHRED']

# Things to be interested in: 
#   1) Error rate of samples relative to expectation.
#   2) Error rate of nucleotides relative to expectation.
#   3) Correspondence between error rate and expectation.
#   4) Distribution of expectations.

from tuba_seq.fastq import get_QC_map
QC_map = get_QC_map()
expected_errors = QC_map[32:][relevant.columns]

samples_relative_to_expectation = stats.groupby(level='Sample').agg(lambda df: ((df['error_rate'] - expected_errors)*df['sum']).sum()/df['sum'].sum())['error_rate']

Log("Mutation rate of samples in excess of PHRED expectations:\n"+samples_relative_to_expectation.to_string(float_format='{:.4%}'.format), True)

def mutation_error_rate(df):
    observed = df.query('Mutation not in @correct')/df.sum()
    return observed - expected_errors/3

full_mut_rate_v_expectation = full_trans.groupby(lambda mut: mut[0]).apply(mutation_error_rate)
L = full_trans.sum()

ave_mut_rate_v_expectation = (full_mut_rate_v_expectation*L).sum(axis=0)/L.sum()

Log("Mutation types in excess of PHRED expectations:\n"+ave_mut_rate_v_expectation.to_string(float_format='{:.4%}'.format), True)

error_rates = 1 - full_trans.apply(lambda col: col.loc[correct].sum()/col.sum())
counts = full_trans.loc[~full_trans.index.isin(correct)].sum()
from scipy.stats.distributions import poisson
interval = 0.95
low_CI, high_CI = poisson.interval(interval, counts)

low_CI *= error_rates/counts
high_CI *= error_rates/counts

from ipy import pdf, plt
ax = plt.gca()
Y = error_rates.values
ax.set_yscale('log')
ax.errorbar(error_rates.index, Y, yerr=(high_CI-Y, Y-low_CI), fmt='.', ecolor='k', elinewidth=2, capsize=0, label='observed')
ax.plot(relevant.columns, expected_errors, color='r', label='PHRED definition')
mu = (L*ave_mut_rate_v_expectation).sum()/L.sum()
ax.text(0.95, 0.95, 'Excess Error Rate: {:.4%}'.format(mu), transform=ax.transAxes, ha='right', va='top')
ax.set(xlabel='Combined PHRED Score', ylabel='Error Rate')
pdf("ErrorRate_v_PHRED")
#plt.savefig('ErrorRate_v_PHRED.pdf', transparent=True, bbox_inches='tight')

for i in range(full_trans.columns.max()):
    full_trans[i] = full_trans[i] if i in full_trans.columns else np.zeros(len(full_trans))
full_trans.sort_index(axis=1).to_csv(args.outfile, index_label=False)
