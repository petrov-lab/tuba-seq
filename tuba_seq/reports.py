import matplotlib
matplotlib.use("Agg")
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

def deconstruct_size_ratios(lSR, min_detectable_contamination, n_components=3 ):
    from sklearn.mixture import GaussianMixture as model
    start = -np.log10(min_detectable_contamination)
    m = model(
        n_components=n_components,
        covariance_type='spherical',
        max_iter=100,               # default: 100
        n_init=1,                   # default: 1
        init_params='random',
        means_init=np.r_[np.zeros(n_components-1), start].reshape(-1, 1)
        )
    return m.fit(lSR.values.reshape(-1, 1))

def contamination(tumor_numbers, alpha=0.5, min_detectable_contamination=0.1):
    from statsmodels.formula.api import ols
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
    mean_BCs = contaminant_BCs.mean()
    def find_contaminants(contaminator, contaminant):
        overlap = contaminant*contaminator
        Size_Ratio = M.loc[overlap, contaminator.name] / M.loc[overlap, contaminant.name]
        OB = overlap.sum()
        EOB = conditional_barcode_freq.loc[contaminator].sum()*contaminant_BCs.loc[contaminant.name]/mean_BCs 
        #print(contaminator.sum(), contaminant.sum(), OB/EOB, poisson.cdf(EOB, OB))
        if OB/EOB > 5:
            lSR = np.log10(Size_Ratio)
            f = plt.figure()
            ax = plt.gca()
            n, bins, patches = ax.hist(lSR, bins=40)
            gmm = deconstruct_size_ratios(lSR, min_detectable_contamination)
            stats = pd.DataFrame(dict(means=gmm.means_.flatten(), weights=gmm.weights_.flatten(), stds=np.sqrt(gmm.covariances_.flatten())))
            bic = gmm.bic(lSR.reshape(-1, 1))
            print(contaminator.name, contaminant.name, 'BIC: {:.1}'.format(bic))
            print(stats.to_string(float_format='{:.2}'.format))
            X = np.linspace(bins[0], bins[-1], 400)
            Y = np.exp(gmm.score_samples(X.reshape(-1, 1)).flatten())
            Y /= Y.max()/n.max()
            ax.plot(X, Y, 'r-') 
            plt.savefig('{:}_{:}.pdf'.format(contaminator.name, contaminant.name), format='pdf', bbox_inches='tight')
        return pd.Series({
            'Fraction_of_Reads'        : M.loc[overlap, contaminant.name].mean() / M.loc[overlap, contaminator.name].mean(),
            'Overlapping_Barcodes'     : overlap.sum(),
            'Expected_Barcode_Overlap' : conditional_barcode_freq.loc[contaminator].sum()*contaminant_BCs.loc[contaminant.name]/mean_BCs,
            'N_contaminant'            : contaminant.sum(), 
            'N_contaminator'           : contaminator.sum()})
        
    df = pd.DataFrame({(name_ator, name_ant):find_contaminants(contaminator, contaminant) for name_ator, contaminator in contaminators.iteritems() for name_ant, contaminant in contaminants.iteritems() if name_ator != name_ant}).T
    df.index.names = ['Contaminator', 'Contaminant']
    return df

    m = N_mice*(N_mice - 1)
    print(df)
    df['p_value'] = df.apply(lambda row: poisson.cdf(row['Expected_Barcode_Overlap'], row['Overlapping_Barcodes']) , axis=1)
    df['FWER'] = df['p_value'].apply(lambda p: min(0.5, p*m))
    regression = ols("Overlaping_Barcodes ~ Expected_Barcode_Overlap", data=df).fit()
    print(regression.summary())
    out = regression.outlier_test(method='sidak', alpha=alpha)
    #out.index.names = ['is contaminated by']
    #out.columns.names = ['Statistic']
    #out['p_value'] = out['unadj_p']*m
    #return out.query('student_resid > 1 and Size_Ratio > 1 and FWER < '+str(alpha))[['FWER', 'Size_Ratio', "Overlapping_Barcodes", 'Expected_Barcode_Overlap']]
    return out.query('Fraction_of_Reads < 1 and FWER < '+str(alpha))[['FWER', 'Fraction_of_Reads', "Overlapping_Barcodes", 'Expected_Barcode_Overlap']]

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
