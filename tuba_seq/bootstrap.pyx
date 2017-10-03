import pandas as pd
from numpy.linalg import norm
import numpy as np
from time import time
from warnings import warn
from pmap import large_iter_pmap as pmap
from itertools import combinations_with_replacement

def possible_boots(n): 
    from scipy.misc import comb
    return comb(n+n - 1, n)

def _rs(df): 
    return df.sample(len(df), replace=True)

def Bonferroni(P, alpha, m):
    return P <= alpha/m

def Sidak(P, alpha, m):
    return P <= 1 - (1 - alpha)**(1/m)

def Holm(P, alpha, m):
    P_ordered = P.sort_values()
    for k, P_k in enumerate(P_ordered):
        if P_k > alpha/(m + 1 - k):
            # R = k
            P_pass = P_k
            break
    return P < P_pass

def Hochberg(P, alpha, m):
    P_ordered = P.osrt_values()
    for k, P_k in enumerate(P_ordered):
        if P_k <= alpha/(m + 1 - k):
            P_failed = P_k
        else:
            break
    return P <= P_failed

class sample(object):
    permissible_uncertainty = 0.21

    def rs(self):
        return _rs(self.df) if self.constrain_level is None else self.gb.apply(_rs)

    def bs_estimate(self, *_args):
        return self.estimator(self.rs(), **self.kwargs)

    def __init__(self, df, estimator, N=4000, time_warning=60, map=pmap, constrain_level=None, **kwargs):
        """sample(df, estimator) -> bootstrap object

Creates bootstrap sampling distribution of df.apply(estimator). 

Input:
    df  : pandas.DataFrame
        Contains data to analysis with estimator
    
    estimator   : function(1-D array-like) -> float
        Operates on each column of df. Can be any estimator e.g. mean, 
        std, median
    
    N (default:2000) : int
        Number of bootstrap samples to create for sampling distribution. In
        general, 1,825 replications are needed to obtain 1% accuracy of 95% CI.
    
    constrain_level (default: None) : Index identifier (int or str)
        Samples an equal number of rows from each category in *constrain_level*.

    **kargs will be passed to 'estimator' function.  

Properties of the Bootstrap Method:
    -Assumes independence of samples (a limitation of nearly all CI estimators)
    -Simple
    -Applicable to any estimator
    -More accurate than CI estimators that assume normality when N is 
    asymptotically large. However, N may need to be very large to converge to 
    true CI, rendering the method infeasible. 
"""
        assert N >= 1, "must have at least one bootstrap sample"
        self.df         = df
        self.estimator  = estimator 
        self.constrain_level = constrain_level
        self.map = map

        self.min_pvalue = 4/(N*self.permissible_uncertainty**2)
        self.kwargs = kwargs
        self.true_estimate = estimator(self.df, **kwargs)
        
        if constrain_level is not None:
            self.gb = df.groupby(level=constrain_level)
        
        if possible_boots(len(df)) < N:
            print("N > len(df)! --> Calculating *entire* bootstrap distribution")
            self.bootstrap_samples = pd.DataFrame([self.estimator(df.loc[list(idxs)], **kwargs) for idxs in combinations_with_replacement(df.index.tolist(), len(df))]) 
        else:
            self.bootstrap_samples = pd.DataFrame(list(map(self.bs_estimate, range(N))))

    def Jackknife(self):
        return np.array([self.estimator(self.df.drop(ix), **self.kwargs) for ix in self.df.index])

    def percentileofscore(self, scores):
        from scipy.stats import percentileofscore
        return (self.bootstrap_samples - scores).apply(percentileofscore, args=(0, 'mean'))
    
    def CI(self, alpha=0.05, bias_corrected=False, for_matplotlib=False):
        """sample.CI(alpha=0.05) -> pandas.df

Estimates Confidence Interval (CI) of any estimator using Bootstrap Method.

Parameters:
    alpha (default:0.05) : float or array-like
        Level of significance or quantiles of sampling distribution to report:
            if float        -> 1 - alpha CI
            if array-like   -> quantiles of sampling distribution
    
    bias_corrected (default:'auto') : bool
        Compensates for skewness in bootstrap distribution. See Efron, B.
        (1987) "Better Bootstrap Confidence Intervals" for details. 
    
        There are a number of sanity-checks and conditions for applying this correction. 

    for_matplotlib (default:False) : bool
        Output can be fed directly into matplotlib errorbar assignment. 
"""
        quantiles = [alpha/2, 1-alpha/2] 
        
        if quantiles[0] < self.min_pvalue:
            warn("There is insufficient sampling for the {:.2%} Confidence Interval.".format(1 - alpha))

        samples = self.bootstrap_samples
        if bias_corrected == 'auto':
            bias_corrected = len(self.df) < len(samples)    # Not useful when the Jackknife takes longer to calculate than the bootstrap distribution (bias scales with len(df)**-0.5, anyways). 
            if bias_corrected:
                print("Correcting bootstrap distribution bias...")
        
        if bias_corrected and not self.constrain_level:    # Only apply the bias correction if bootstrapping was on individual samples 

            Jackknife = self.Jackknife()
            Theta = Jackknife - Jackknife.mean()
            a = (Theta**3).sum(axis=0) / ( 6*((Theta**2).sum(axis=0)**1.5) )
            z_true = norm.ppf(self.percentileofscore(self.true_estimate)*1e-2)[:, np.newaxis]
            z_q = z_true + norm.ppf(quantiles)
            adjusted_z = z_true + z_q/(1 - a[:, np.newaxis]*z_q)
            adjusted_quantiles = norm.cdf(adjusted_z) 
            
            assert (adjusted_quantiles[:,0] <= 0.5).all() and (adjusted_quantiles[:,1] >= 0.5).all(), "Bias Correction gave absurd adjusted_quantiles."
            CIs = pd.DataFrame([S.quantile(q=quantiles).values for quantiles, (_col, S) in zip(adjusted_quantiles, samples.iteritems()) ], columns=['low', 'high'], index=samples.columns)
        else:
            CIs = samples.quantile(q=quantiles).T
            CIs.columns = 'low', 'high'

        above_true = CIs['low'] > self.true_estimate
        if above_true.any():
            warn("Low CI above True Estimate.")
            CIs.loc[above_true, 'low'] = np.nan
        
        below_true = CIs['high'] < self.true_estimate
        if below_true.any():
            warn("High CI below True Estimate.")
            CIs.loc[below_true, 'high'] = np.nan
        
        if for_matplotlib:
            return np.vstack((self.true_estimate - CIs['low'], CIs['high'] - self.true_estimate)) 
        
        CIs.insert(len(alpha)/2 if hasattr(alpha, '__len__') else 1, 'true', self.true_estimate)
        CIs.columns.name = 'Estimates' 
        return CIs 

    def print_CI(self, alpha=0.05, formatter='.2f', **kargs):
        CIs = self.CI(alpha, **kargs)
        return '\n'.join(["{name:<10}: {true: FORMAT} {0:.0%} CI [{low: FORMAT}, {high: FORMAT}]".replace('FORMAT', formatter).format(
                        1 - alpha, **locals()) for name, (low, true, high) in CIs.iterrows()])

    def pscores(self, null_hypothesis=0, two_sided=True):
        """null_hypothesis [def: 0] : value to reject.
        test [def: 'two-sided']  : can also be 'above' or 'below' null_hypothesis.
        bonferroni_correction [def: 1] : Number of multiple-hypotheses
"""
        pscores = self.percentileofscore(null_hypothesis)*1e-2
        return pd.DataFrame([pscores, 1 - pscores]).min()*2 if two_sided else pscores

    def pstars(self, alphas=np.array([0.05, 0.01, 0.001, 0.0001]), FWER_method=Bonferroni, **kwargs):
        """null_hypothesis [def: 0] : value to reject.
        test [def: 'two-sided']  : can also be 'above' or 'below' null_hypothesis.
        bonferroni_correction [def: 1] : Number of multiple-hypotheses
"""
        m = kwargs.pop('m', len(self.true_estimate))

        P = self.pscores(**kwargs)
        H = pd.concat([FWER_method(P, alpha, m) if FWER_method is not None else (P <= alpha) for alpha in alphas], axis=1)
         
        return H.sum(axis=1).apply(lambda n: n*'*')

def describe(df, estimator, **kargs):
    """ describe(df, estimator_string) -> Readable CI for statstic_str
"""
    obj = sample(df, estimator, **kargs) 
    return "Bootstrap estimates of {:}:".format(estimator.__name__ if hasattr(estimator, '__name__') else 'estimator') + '\n' + obj.print_CI()

def seaborn_hacker(S, axis=None):
    if axis is None:
        return S[1]
    if axis == 1:
        return S[:, 1]
    if axis == 0:
        return S[1, :]
   
def rsquared(df, levels=0, agg_level=None):
    """"""
    from statsmodels.api import OLS
    expanded = df.reset_index(levels)
    expanded.index = df.index
    expanded = expanded.dropna()
    condensed = expanded.groupby(level=agg_level).agg(np.mean) if agg_level else expanded
    Y = condensed.pop(condensed.columns[-1])
    result = OLS(Y, condensed).fit()
    return result.rsquared

if __name__ == '__main__':
    np.random.seed(seed=42)
    N = 1000
    x = pd.DataFrame(np.random.randn(N, 3) + np.arange(3), columns=['dogs', 'cats', 'monkeys'])
    print(describe(x, np.mean))
    print(describe(x, np.std ))
