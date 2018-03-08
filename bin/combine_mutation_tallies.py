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

correct = {'A2A', 'C2C', 'G2G', 'T2T'}

data = pd.concat({f.stem:pd.read_csv(f, index=[0]) for f in args.input_dir.iterdir() if f.suffixes[0] == '.csv'}, names=['Sample'])
full_trans = data.groubpy(level='mutation').sum()
full_trans.to_csv(args.outfile)

# Now the report
gb = data.groupby(level=['Sample'])

def error_rate(col):
    return 1 - col.loc[correct].sum()/col.sum()

stats = gb.agg([error_rate, len])

# Things to be interested in: 
#   1) Error rate of samples relative to expectation.
#   2) Error rate of nucleotides relative to expectation.
#   3) Correspondence between error rate and expectation.
#   4) Distribution of expectations.

from tuba_seq.fastq import QC_map
expected_errors = QC_map[32:][:data.columns.max()]

def error_relative_to_expectation(df):
    L = df.loc['len']
    N = L.sum()
    observed = (df.loc['error_rate']*L).sum()/N
    expected = (expected_errors[df.columns]*L)/N
    return observed - expected

samples_relative_to_expectation = stats.stack().stack().groupby(level='Sample').agg(error_relative_to_expectation)

Log("Mutation rate of samples in excess of PHRED expectations:\n"+samples_relative_to_expectation.to_string(float_format='{:.4%}'.format), True)

def mutation_error_rate(col):
    gb = col.groupby(lambda mut: mut[0])
    rates = gb.transform(lambda f: f/f.sum())
    rates -= expected_errors[col.name]/3 
    rates.loc[correct] = np.nan
    return rates

all_samples = gb.sum()
full_mut_rate_v_expectation = all_samples.apply(mutation_error_rate)
L = all_samples.sum()
ave_mut_rate_v_expectation = (full_mut_rate_v_expectation*L).sum(axis=0)/L.sum()

Log("Mutation types in excess of PHRED expectations:\n"+ave_mut_rate_v_expectation.to_string(float_format='{:.4%}'.format), True)

error_rates = all_samples.apply(error_rate)
counts = all_samples.loc[~all_samples.index.isin(correct)].sum()

from scipy.stats.distribution import poisson
interval = 0.95
low_CI, high_CI = poisson.interval(interval, counts)*error_rates/counts

import matplotlib.pyplot as plt
ax = plt.gca()
Y = error_rates.values
ax.set_yscale('log')
ax.errorbar(error_rates.index, Y, yerr=(high_CI-Y, Y-low_CI), fmt='.', ecolor='k', elinewidth=2, capsize=0, label='observed')
ax.plot(np.arange(len(expected_errors)), expected_errors, color='r', label='PHRED definition')
mu = (error_rates*counts).sum()/counts
ax.text(0.95, 0.95, 'Error Rate: {:.4%}'.format(mu), transform=ax.transAxes, ha='right', va='top')
ax.set(xlabel='Combined PHRED Score', ylabel='Error Rate')
plt.savefig('ErrorRate_v_PHRED.pdf', transparent=True, bbox_inches='tight')

