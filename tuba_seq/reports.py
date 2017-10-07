import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tuba_seq.graphs import text_color_legend
from scipy import stats 

def contamination(tumor_numbers, alpha=0.05):
    from statsmodels.formula.api import ols
    
    M = tumor_numbers.unstack(level='Mouse')
    I = M.notnull()
    N_mice = len(I.columns)
    m = N_mice*(N_mice - 1)/2

    def find_contaminants(S_m):
        N = S_m.sum()
        I_other = I.drop(S_m.name, axis=1)
        Overlap = I_other.mul(S_m, axis=0)
        N_overlap = Overlap.sum()
        N_both = I_other.sum() + N - N_overlap
        df = pd.DataFrame(dict(N_both=N_both, enrichment=N_overlap/N_both))
        regression = ols("enrichment ~ N_both", data=df).fit()
        out = regression.outlier_test()
        out.index.names = ['is contaminated by']
        out.columns.names = ['Statistic']
        out['p_value'] = out['unadj_p']*m
        out['Size_Ratio'] = Overlap.apply(lambda Overlap_i: M.loc[Overlap_i.values, S_m.name].mean() / M.loc[Overlap_i.values, Overlap_i.name].mean(), axis=0)
        return out.query('student_resid > 1 and Size_Ratio > 1 and p_value < '+str(alpha))[['p_value', 'Size_Ratio']]
    return pd.concat({mouse:find_contaminants(S_m) for mouse, S_m in I.iteritems()}, names=['Mouse'])

def barcode_diversity(tumor_numbers, plot=True):
    sns.set_style('whitegrid')
        # Frequency profile
    barcodes = tumor_numbers.index.get_level_values('barcode')
    random_barcode_length = len(barcodes.values[0])
    PWM = pd.concat({i:barcodes.str.get(i).value_counts()/len(barcodes) for i in range(random_barcode_length)})
    PWM.index.names = ['Location', 'Nucleotide']
    PWM.name = 'Frequency'
    if plot:
        f, (ax1, ax2, ax3) = plt.subplots(3, figsize=(6, 12))
        base = np.zeros(random_barcode_length)
        for nuc in 'ACGT':
            X = PWM.loc[:, nuc]
            ax1.bar(X.index.values, X, bottom=base, label=nuc)
            base += X.values
        yticks = [0, 0.25, 0.5, 0.75, 1]
        ax1.set(ylabel='Frequency', xlabel='Position', ylim=[0, 1], yticks=yticks, yticklabels=list(map('{:.0%}'.format, yticks)))
        text_color_legend(ax1, title='Nucleotide:') 
        sns.despine(left=True, right=True, top=False)

    def sgID_statistics(S):
        M = S.unstack(level='Mouse').notnull()
        observed = len(M)
        barcodes_per_mouse = M.sum()
        N = len(barcodes_per_mouse)
        mu_geq1 = M.sum(axis=1).mean()
        from scipy.optimize import brentq
        def P_0(p):
            return (1 - p)**N
        p = brentq(lambda p: mu_geq1 - p*N/P_0(p), 0, mu_geq1/N, full_output=False)
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
