#!/usr/bin/env python3
import argparse, os
import pandas as pd
import numpy as np
from tuba_seq.shared import logPrint
from tuba_seq.reports import plt
import seaborn as sns

tuba_seq_dir = os.path.dirname(__file__)

############################## Input #########################################

parser = argparse.ArgumentParser(   description="Corrects consolidated & annotated clustering output for GC-bias and transforms abundances to absolute cell number, and applies final filter. Also, can provide a final report on the diversity of the random barcodes.", 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

IO_group = parser.add_argument_group('IO', 'Input/Output Optional Arguments')
IO_group.add_argument('-i', '--input_file', default='combined.csv.gz', help='CSV file with the consolidated samples and their sgID-target annotations.')
IO_group.add_argument('--full_out_file', type=str, default='complete_analysis.csv', help='Comprehensive description of every un-dropped barcode cluster.')
IO_group.add_argument('--simple_out_file', type=str, default='tumor_sizes.csv', help='Only tumors that pass the filter and their corrected sizes.') 

IO_group.add_argument('-b', '--bartender', action='store_true', help='Process Bartender (https://github.com/LaoZZZZZ/bartender-1.1) clustering output.')
IO_group.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
IO_group.add_argument('-s', '--spike_barcodes', nargs='+', default='infer', help='Barcodes present in the spiked-in benchmark barcodes')
IO_group.add_argument('--report', action='store_true', help='Graph the best-fitting model of GC bias, and provide a report on Barcode Diversity & potential contaminations.')

OP_group = parser.add_argument_group('OP', 'Optional arguments affecting operation')
OP_group.add_argument('--final_filter', default='Concentration >= 1 and Cells >= 500', help='Final barcode filter')
OP_group.add_argument('-p', '--parallel', action='store_true', help='Parallelize operation.')
OP_group.add_argument('--spike_cells', type=float, default=500000, help='Number of cells in the spiked-in barcodes.')
OP_group.add_argument('-c', '--cell_metric', type=str, default='n0 + n1', help='Measure of cell abundance to use (can be any valid pandas.DataFrame.eval argument).')
OP_group.add_argument('--proportion', type=float, default=0.01, 
    help='Proportion of values to trim for trimmed mean (this was zero in the original paper, but the default value is now recommended).')
OP_group.add_argument("--find_order", action="store_true", help="Identifies the best-fitting polynomial function.")
OP_group.add_argument('-o', '--order', type=int, default=4, help='Polynomial order to use, if --find_order is not invoked.')
OP_group.add_argument('--max_order', type=int, default=6, help='Maximum Order of linear-fit to consider.')
OP_group.add_argument('-d', '--drop', type=str, nargs='+', default=('Unknown',), help='sgRNAs to drop before fitting data for GC-bias correction.')
OP_group.add_argument('-l', '--linear', action='store_true', help="Will not log-transform data before fitting GC-content (not recommended)")

args = parser.parse_args()
Log = logPrint(args)
if args.parallel:
    from tuba_seq.pmap import pmap as map 

combined = pd.read_csv(args.input_file, index_col=[0, 1])
clean = combined.query("target not in @args.drop")

########################## GC-Correction ######################################

def trimmed_mean(S):
    from scipy.stats import trim_mean
    return trim_mean(S, args.proportion) if args.proportion > 0 else np.mean(S)

clean.insert(0, 'GCs', clean['ID'].str.count("G") + clean['ID'].str.count("C") + clean['barcode'].str.count("G") + clean['barcode'].str.count("C"))

barcode_length = clean.iloc[0][['sgRNA', 'barcode']].str.len().sum()

clean.set_index(['barcode', 'GCs'], append=True, inplace=True)
untransformed = clean.eval(args.cell_metric)
Y = untransformed if args.linear else np.log(untransformed)
residual = Y.groupby(level=['Mouse', 'target']).transform(lambda S: S - S.mean())
gb = residual.groupby(level='GCs')

data = gb.agg(trimmed_mean)

def model_fit(order, x=data.index.values, Y=data.values, weights=gb.count()):
    import statsmodels.formula.api as sm
    X = pd.DataFrame({'x^{:}'.format(i):x**i for i in range(order+1)}, index=x)
    return sm.WLS(Y, X, weights=weights).fit()

if args.find_order:
    orders = np.arange(2, args.max_order+1)
    models = map(model_fit, orders)
    R2 = pd.Series([model.rsquared_adj for model in models], index=orders)
    args.order = R2.argmax()
    Log.line_break()
    Log("Comparison of models of GC-bias")
    Log.line_break()
    Log('Order | Adj. R2 | To be used')
    for order, r2_i in R2.iteritems():
        Log("{:}      {:.2%}      {:}".format(order, r2_i, '*' if order == args.order else ''))

model = model_fit(args.order)
Log(model.summary(), True)

if args.report:
    df = residual.reset_index()
    X_lim = df['GCs'].quantile(q=[0.0005, 0.9995])
    x = np.linspace(*X_lim.values, num=200)
    X = pd.DataFrame({'x^{:}'.format(i):x**i for i in range(args.order+1)}, index=x)
    plot_data = df.query("GCs >= {:} and GCs <= {:}".format(*X_lim.values))
    ax = sns.pointplot(x='GCs', y=0, estimator=lambda x: np.exp(trimmed_mean(x)) if not args.linear else trimmed_mean(x), data=plot_data, join=False)
    P = model.predict(X.values) if args.linear else np.exp(model.predict(X.values))
    ax.plot(x - x[0], P)
    ax.text(0.95, 0.95, 'Adjusted $R^2$ = {:.0%}'.format(model.rsquared_adj), transform=ax.transAxes, ha='right', va='top')
    ax.set(xlabel='GC Content', ylabel='Marginal Effect on Tumor Size')
    
        # Convert x-tick labels to % GC content (rather than GC tallies) 
    xrange_frac = np.arange(0, 1, 0.1)
    xrange_GCs = xrange_frac*barcode_length
    keep = (xrange_GCs > X_lim.iloc[0])&(xrange_GCs < X_lim.iloc[1])
    ax.xaxis.set_ticks(xrange_GCs[keep] - x[0])
    ax.xaxis.set_ticklabels(list(map(str, xrange_frac[keep])))
    plt.savefig("Quality_of_GC_fit", transparent=True, bbox_inches='tight')

    

predictor = pd.Series(model.predict(), index=data.index)

predictions = predictor.loc[residual.index.get_level_values('GCs')].values
clean.insert(0, 'GC_corrected', Y - predictions if args.linear else np.exp(Y - predictions))

################# Calculate Absolute Cell Number & final Statistics ###########

if args.spike_barcodes == 'infer':
    min_necessary_spike_fraction = 0.01    #fraction of largest barcode to qualify as legitimate
    spike_candidates = clean.loc[(slice(None), 'Spike'), 'GC_corrected'].groupby(level='barcode').sum()
    candidate_frac = spike_candidates/spike_candidates.sum()
    candidate_frac.name = 'fractional abundance'
    spike_BCs = candidate_frac.loc[candidate_frac > min_necessary_spike_fraction]
    Log("Identified the following spike-in barcodes as legit:", True)
    Log(spike_BCs, True)
    spike_barcodes = spike_BCs.index.values
else:
    spike_barcodes = args.spike_barcodes

true_spikes = clean.query('target == "Spike" and barcode in @spike_barcodes')

scalings = true_spikes['GC_corrected'].groupby(level='Mouse').mean()

def spike_normalization(df):
    return df*args.spike_cells/scalings.loc[df.name]

final = clean.join(pd.DataFrame(
       {'Cells': clean['GC_corrected'].groupby(level='Mouse').transform(spike_normalization),
        'Concentration':clean.eval('n0 / (abundance - n1 - n0)')}))
final.insert(0, 'pass_filter', final.eval(args.final_filter))
final.insert(0, 'false_positive', final.eval('target == "Spike" and pass_filter and not barcode in @spike_barcodes'))

# Output 

false_positives = final['false_positive'].sum()
FDR = false_positives/(false_positives + len(true_spikes))

Log("""Estimated FDR: {:.1%} (Based on the number of unexpected spike-in barcodes that passed the final filter.)""".format(FDR, args.final_filter), True)

final.to_csv(args.full_out_file+'.gz', compression='gzip')
tumor_df = final.query('pass_filter and target != "Spike"').reset_index()[['Mouse', 'target', 'barcode', 'Cells']]
tumor_df.to_csv(args.simple_out_file+'.gz', index=False, compression='gzip')
tumors = tumor_df.set_index(["Mouse", 'target', 'barcode'])["Cells"]

if args.report:
    from tuba_seq.reports import barcode_diversity, contamination
    Log("Graphing on the Diversity of the Random Barcodes...")
    barcode_diversity(tumors)
    contaminants = contamination(tumors) 
    if len(contaminants) == 0:
        Log("Found no evidence of cross-contamination of samples.", True)
    else:
        Log("""
There was an abnormally-high degree of overlapping barcodes between specific 
samples in this library suggesting cross-contamination of the following samples:
"""+contaminants.reset_index().to_string(index=False)+"""
Directionality of contamination is inferred by the `Size Ratio` of overlapping
barcodes, i.e. the contamination volume should be smaller than its original 
sample. This may be an unreliable assumption.""", True)

