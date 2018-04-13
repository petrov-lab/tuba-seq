import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from tuba_seq.graphs import text_color_legend
from tuba_seq.tools import LN_mean, inerts
from scipy import stats 

def identify_outliers_by_target_profile(normalized_tumors, metric=LN_mean, alpha=0.05, inerts=inerts, sample_level='Mouse'):
    from tuba_seq.bootstrap import sample
    from scipy.stats import combine_pvalues

    expected = normalized_tumors.groupby(level='target').agg(metric)

    def mouse_specific(t):
        return t.groupby(level=[sample_level, 'target']).agg(metric).groupby(level=sample_level).transform(lambda S: S - expected.values)
        
    m = len(normalized_tumors.groupby(level=['target', sample_level]))
    
    bs = sample(normalized_tumors, mouse_specific, min_pvalue=alpha/m)

    pscores = bs.pscores(null_hypothesis=0, two_sided=True)
    pvals = pscores.groupby(level=sample_level).agg(lambda S: combine_pvalues(S)[1])*m
    outliers = bs.true_estimate.unstack(level='target').loc[pvals < alpha]
    if len(outliers) == 0:
        print("Did not find outliers")
        return None

    outliers.insert(0, 'p-value', pvals.loc[pvals < alpha])

    untransformed_outlier_bootstraps = bs[pvals < alpha].T.groupby(level=sample_level).transform(lambda S: S + expected.values).T
        # Since we calculated all these bootstraps, I'm going to hack the bootstrap object to see if any of the active
        # sgRNAs exhibit increased tumor size.
    bs.bootstrap_samples = untransformed_outlier_bootstraps
    
    most_active_sgRNA_pscore = bs.pscores(null_hypothesis=1, two_sided=False).groupby(level=sample_level).agg(
                                    lambda S: S.loc[-S.index.get_level_values("target").isin(inerts)].min())
    
    active_m = (-expected.index.isin(inerts)).sum()
    outliers.insert(1, 'Cas9 Activity?', most_active_sgRNA_pscore.apply(lambda p: 'yes' if p < alpha/active_m else 'no'))
    return outliers

inverse_dict={np.sqrt:np.square, np.log:np.exp}
def deconstruct_size_ratios(SR, min_detectable_contamination=0.05, n_components=3, transform=np.sqrt, std_threshold=0.25):
    from sklearn.mixture import GaussianMixture as model
    inv = inverse_dict[transform]
    tSR = transform(SR)
    m = model(
        n_components=n_components,
        covariance_type='spherical',
        max_iter=1000,                # default: 100
        n_init=50,                   # default: 1
        init_params='random')
    gmm = m.fit(tSR.values.reshape(-1, 1))
    stats = pd.DataFrame(dict(mean=gmm.means_.flatten(), weight=gmm.weights_.flatten())).apply(inv)
    stats['std'] = np.sqrt(gmm.covariances_.flatten())
    contaminants = stats.query('mean > 1 and mean < 1/@min_detectable_contamination and weight > @min_detectable_contamination and std <= @std_threshold')
    bic = gmm.bic(tSR.values.reshape(-1, 1))
    return dict(isContamination=len(contaminants) == 1 and bic < 0,
                stats=stats, 
                Contaminating_Fraction = 1/(1 + contaminants['mean'].values[0]) if len(contaminants) == 1 else np.nan, 
                GMM=gmm,
                bic=bic, 
                transform=transform)

def contamination(tumor_numbers, alpha=0.01, min_detectable_contamination=0.05, map=map, graph_contaminations=True, min_size_ratio=0.75, graph_columns=4):
    from scipy.stats.distributions import poisson
    M = tumor_numbers.unstack(level='Sample')
    min_tumor_size = tumor_numbers.min()
    contaminators = M.apply(lambda S: S > min_tumor_size/min_detectable_contamination)
    contaminants = M.notnull()
    N_mice = len(M.columns)
    if N_mice < 3:
        print("Can't find contaminants with less than 3 mice.")
        return pd.DataFrame()
    conditional_barcode_freq = (contaminants.sum(axis=1) - 1)/(N_mice - 1)
    contaminant_BCs = contaminants.sum()
    mean_contaminants = contaminant_BCs.mean()
    contamination_propensity = (conditional_barcode_freq*contaminators.T).sum(axis=1)
    df = pd.DataFrame({
        'Overlapping Barcodes'     : contaminators.T.dot(contaminants.astype(int)).stack(),
        'Expected Barcode Overlap' : contamination_propensity.apply(lambda x: x*contaminant_BCs/mean_contaminants).stack()})

    df.index.names = ['Contaminator', 'Contaminant']
    df = df.query('Contaminator != Contaminant')
    #from statsmodels.formula.api import ols
    #regression = ols("Overlapping_Barcodes ~ Expected_Barcode_Overlap + N_contaminant + N_contaminator", data=df).fit()
    #print(regression.summary())
    #df['Empirical_Expected_Barcode_Overlap'] = regression.predict()
    df.insert(0, 'P-value', df.apply(lambda row: poisson.cdf(row['Expected Barcode Overlap'], row['Overlapping Barcodes']) , axis=1))
    m = N_mice*(N_mice - 1)
    df.insert(0, 'FWER', df['P-value'].apply(lambda p: min(0.5, p*m)))
    contaminations = df.query('FWER < @alpha')
    print('{:.1%} ({:} sample pairs) have anomalously-high fractions of overlapping-barcodes...'.format(
          len(contaminations)/len(df), len(contaminations)))
    
    ix = contaminations.index
    overlaps = (contaminants[ix.get_level_values('Contaminant')].values&contaminators[ix.get_level_values("Contaminator")].values)
    
    contaminations.insert(0, "Size Ratios", [(M.loc[overlap, contaminator]/M.loc[overlap, contaminant]).loc[lambda x: x >= min_size_ratio] for (contaminator, contaminant), overlap in zip(ix, overlaps.T)])
    from functools import partial
    f = partial(deconstruct_size_ratios, min_detectable_contamination=min_detectable_contamination)
    
    
    
    outputs = pd.DataFrame(map(f, contaminations['Size Ratios']), index=ix)
    transform = outputs.pop('transform').values[0] 
    final = pd.concat([contaminations, outputs], axis=1)
    final = final.loc[final.pop('isContamination')]
    
    if graph_contaminations: 
        Max = final['Size Ratios'].apply(max).max()
        bins = np.geomspace(min_size_ratio, Max, 50)
        X = np.geomspace(min_size_ratio, Max, 400)
        with sns.axes_style('white'):
            L = len(contaminations)
            graph_rows = int(np.ceil(L/graph_columns))
            f, axs = plt.subplots(graph_rows, graph_columns, sharex=True, sharey=True, figsize=(graph_rows*4, graph_columns*4))
            axs = axs.flatten()
            for ax, (ix, row) in zip(axs.flatten(), final.iterrows()):
                n, __, patches = ax.hist(row['Size Ratios'], bins=bins)
                Y = np.exp(row['GMM'].score_samples(transform(X).reshape(-1, 1)).flatten())
                Y /= Y.max()/n.max()
                ax.plot(X, Y, 'r-')
                ax.set_title('{:} / {:}'.format(*ix))
            ax.set_xscale('log')
            ax.set(ylabel='# of Overlapping Barcodes', xlabel='Ratio of barcode abundances')
            plt.savefig('contaminations_GMM.pdf', format='pdf', bbox_inches='tight')
    
    return final.drop(['stats', 'Size Ratios', 'GMM'], axis=1)

def barcode_diversity(tumor_numbers, plot=True):
    sns.set_style('whitegrid')
        # Frequency profile
    barcodes = tumor_numbers.index.get_level_values('barcode')
    random_barcode_length = len(barcodes.values[0])
    PWM = pd.DataFrame({i:barcodes.str.get(i).value_counts() for i in range(random_barcode_length)}).fillna(0)/len(barcodes)
    if plot:
        f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6, 12))
        base = np.zeros(random_barcode_length)
        for nuc in 'ACGT':
            X = PWM.loc[nuc]
            ax1.bar(X.index.values, X, bottom=base, label=nuc)
            base += X.values
        yticks = [0, 0.25, 0.5, 0.75, 1]
        ax1.set(ylabel='Frequency', xlabel='Position', ylim=[0, 1], yticks=yticks, yticklabels=list(map('{:.0%}'.format, yticks)))
        text_color_legend(ax1, title='Nucleotide:') 
        sns.despine(left=True, right=True, top=False)

    def sgID_statistics(S):
        M = S.unstack(level='Sample').notnull()
        observed = len(M)
        barcodes_per_mouse = M.sum()
        N = len(barcodes_per_mouse)
        mu_geq1 = M.sum(axis=1).mean()
        from scipy.optimize import brentq
        def P_0(p):
            return (1 - p)**N
        p = brentq(lambda p: mu_geq1 - p*N/P_0(p), 0, mu_geq1/N - np.finfo(np.double).eps, full_output=False)
        unobserved = observed*P_0(p)/(1 - P_0(p))

        def collision_probability(N):
            P_0, P_1 = stats.distributions.poisson.pmf([0, 1], N/(observed+unobserved))
            return (1 - P_0 - P_1)/(1 - P_0)
        
        P_collision_by_mouse = barcodes_per_mouse.apply(collision_probability)
        return {'Estimated Unobserved Barcodes' : unobserved,
                'Observed Barcodes' : observed,
                'Probability of Barcode Collision' : P_collision_by_mouse.dot(barcodes_per_mouse)/barcodes_per_mouse.sum()}

    sgID_info = tumor_numbers.groupby(level='target').apply(sgID_statistics).unstack()
    if plot:
        X = np.arange(len(sgID_info))
        bars = ax2.bar(X, sgID_info['Observed Barcodes'], label='Observed')
        ax2.bar(X, sgID_info['Estimated Unobserved Barcodes'], label='Estimated Unobserved',
                bottom=sgID_info['Observed Barcodes'], color=bars[0].get_facecolor(), alpha=0.5)
        ax2.set_ylabel('Number of Barcodes')
        text_color_legend(ax2, title='Barcodes:', bbox_to_anchor=(0.9, 1)) 
        
        ax3.bar(X, sgID_info['Probability of Barcode Collision'])
        max_ytick = int(np.ceil(sgID_info['Probability of Barcode Collision'].max()*100))
        yticks = np.arange(max_ytick+1)/100
        ax3.set(ylabel='Probability of Barcode Collision', yticks=yticks, ylim=[0, max_ytick/100], yticklabels=list(map('{:.0%}'.format, yticks)))
        for ax in (ax2, ax3):
            ax.set(xticks=X, xlabel='Target')
            sns.despine(left=True, right=True, top=False, ax=ax)
            ax.xaxis.set_ticklabels(sgID_info.index.values, rotation=90, style='italic')
            ax.xaxis.grid(False)

        plt.savefig('barcode_diversity_report.pdf', format='pdf', bbox_inches='tight')
    return sgID_info
