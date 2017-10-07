from numpy cimport ndarray 
from tuba_seq.pmap import pmap as map 
import numpy as np
import pandas as pd
import io, sys, collections

inerts = ['Neo1', 'Neo2', 'Neo3', 'NT1', 'NT3']
percentile_tiers = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99])

def LN_mean(data):
    cdef:
        ndarray[dtype=double, ndim=1] x = np.array(data)
        double L = len(x)
        ndarray[dtype=double, ndim=1] One = np.ones_like(x)
        ndarray[dtype=double, ndim=1] LN_x = np.log(x)
    X = LN_x.dot(One)/L
    X2 = LN_x.dot(LN_x)/L - X*X
    return np.exp(X + 0.5*X2)

def percentiles(S, tiers=percentile_tiers):
    out = S.quantile(q=tiers) 
    final = pd.Series(out.values, index=out.index.map(lambda s: '{:2g}'.format(s*100)))
    final.index.names = ['Percentile']
    return final

def best_power_law_fit(S, resolution=400, x_quantile_range=[0.1, 0.99], sigma_threshold=0.05):
    import powerlaw as pl
    X_lims = S.quantile(q=x_quantile_range)
    X_mins = np.linspace(*X_lims, num=resolution)
    def f(xmin):
        fit = pl.Fit(S, xmin=xmin, sigma_threshold=sigma_threshold, fit_method='Likelihood')
        log_likelihood = fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)[0] 
        return fit, log_likelihood
    Fits = pd.DataFrame(map(f, X_mins), columns=['fit', 'LL'])
    return Fits.loc[Fits['LL'].argmax()]

def best_parametric_distribution(S, consider=['norm', 'expon', 'logistic', 'gumbel', 'lognorm', 'loglogistic']):
    from scipy.stats import anderson
    fits = {dist:anderson(np.array(np.log(S) if dist[:3] == 'log' else S), dist=dist.replace('log', '')).statistic for dist in consider}
    return min(fits, key=fits.get)

class CapturingDict(collections.Counter):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = io.StringIO()
        return self
    
    def __exit__(self, *args):
        self.update(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

def distribution_summarize(S, PL_pval_threshold=0.05, **kargs):
    fit, LL = best_power_law_fit(S, **kargs)
    #with CapturingDict() as output:
    #    ll = fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)[0] # This is a way of suppressing the copious output of this function
    #for k, v in output.items(): print(v, ":", k)
    
    pvalue = 1/(1+np.exp(LL))
    return pd.Series({
        'best_fit':'LN+PL' if pvalue < PL_pval_threshold else 'LN',
        'log_likelihood': LL,
        'alpha':fit.power_law.alpha, 
        'sigma':fit.power_law.sigma, 
        'x_min':fit.power_law.xmin, 
        'P-value':pvalue,
        'LN_Mean':LN_mean(S),
        'log_likelihood_of_truncated_exponential_fit':fit.distribution_compare('power_law', 'truncated_power_law', normalized_ratio=True)[0],
        'Fraction_of_tumors_in_PL':(S > fit.xmin).mean()})

from sklearn.decomposition import PCA
def PC1_weights(cohort, model=PCA(n_components=1)):
    A = cohort.unstack("target")
    model.fit(A)
    weights = model.score_samples(A)
    return pd.concat({rna:pd.Series(weights/weights.sum(), index=A.index) for rna in A.columns}, names=['target'])

def PC1_explained_variance(S, pca=PCA()):
    M = S.unstack("target")
    pca.fit(M)
    return pca.explained_variance_ratio_[0]

def series_ix(S, **kargs):
    return S.loc[(kargs.get(ix, slice(None)) for ix in S.index.names)] 

def inert_normalize(S, estimator='mean', inerts=inerts):
    N = S.loc[S.index.get_level_values('target').isin(inerts)].agg(estimator)
    return S/N

def identify_outliers_by_target_profile(normalized_tumors, metric=LN_mean, method='Mahalanobis', alpha=0.05, inerts=inerts):
    from tuba_seq.bootstrap import sample
    from scipy.stats import combine_pvalues

    expected = normalized_tumors.groupby(level='target').agg(metric)

    def mouse_specific(t):
        return t.groupby(level=['Mouse', 'target']).agg(metric).groupby(level='Mouse').transform(lambda S: S - expected.values)
        
    m = len(normalized_tumors.groupby(level=['target', 'Mouse']))
    
    bs = sample(normalized_tumors, mouse_specific, min_pvalue=alpha/m)

    pscores = bs.pscores(null_hypothesis=0, two_sided=True)
    pvals = pscores.groupby(level='Mouse').agg(lambda S: combine_pvalues(S)[1])*m
    outliers = bs.true_estimate.unstack(level='target').loc[pvals < alpha]
    if len(outliers) == 0:
        print("Did not find outliers")
        return None

    outliers.insert(0, 'p-value', pvals.loc[pvals < alpha])

    untransformed_outlier_bootstraps = bs[pvals < alpha].T.groupby(level='Mouse').transform(lambda S: S + expected.values).T
        # Since we calculated all these bootstraps, I'm going to hack the bootstrap object to see if any of the active
        # sgRNAs exhibit increased tumor size.
    bs.bootstrap_samples = untransformed_outlier_bootstraps
    most_active_sgRNA_pscore = bs.pscores(null_hypothesis=1, two_sided=False).groupby(level='Mouse').agg(
                                    lambda S: S.loc[-S.index.get_level_values("target").isin(inerts)].min())
    
    outliers.insert(0, 'Cas9 Activity?', most_active_sgRNA_pscore.replace({True:'yes', False:'no'}))

    return outliers
