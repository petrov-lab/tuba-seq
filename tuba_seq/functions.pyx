from numpy cimport ndarray 
import numpy as np
import pandas as pd
import io, sys, collections

inerts = ['Neo1', 'Neo2', 'Neo3', 'NT1', 'NT3']
percentile_tiers = np.array([50, 60, 70, 80, 90, 95, 99])

def LN_mean(data):
    """MLE of mean of data, presuming a LogNormal Distribution."""
    cdef:
        ndarray[dtype=double, ndim=1] x = np.array(data)
        double L = len(x)
        ndarray[dtype=double, ndim=1] One = np.ones_like(x)
        ndarray[dtype=double, ndim=1] LN_x = np.log(x)
    X = LN_x.dot(One)/L
    X2 = LN_x.dot(LN_x)/L - X*X
    return np.exp(X + 0.5*X2)

def _xformed_LN_mean(ln_S):
    """LN Mean of vector that is already log-transformed."""
    cdef:
        ndarray[dtype=double, ndim=1] LN_x = ln_S.values
        double L = len(LN_x)
        ndarray[dtype=double, ndim=1] One = np.ones_like(LN_x)
    X = LN_x.dot(One)/L
    X2 = LN_x.dot(LN_x)/L - X*X
    return np.exp(X + 0.5*X2)

def _all_targets_LN_mean(ln_S):
    """LN Mean for all sgRNA targets from Log-transformed vector of tumor sizes."""
    return ln_S.groupby(level='target').agg(_xformed_LN_mean)

def LN_mean_P_values(S, inerts=inerts, min_pvalue=0.0001):
    """Hypothesis test of increased growth relative to sgInerts.

    One-sided bootstrapping hypothesis test that the LN Mean of each sgRNA's 
    tumor size distribution is different from the sgInert distribution using 
    the actual sgInert size distributions as the null sampling distribution. 
    This approach is slightly different than the statistical test used in 
    Rogers et al (2017). See tuba_seq.tools.LN_Mean_Summary for more details. 

    Parameter:
    ----------
    inerts : sgRNAs that are not expected to alter tumor sizes to be used as
             the null sampling distribution (default: ['Neo1-3', 'NT1', 'NT3'])

    min_pvalue : Maximum statistical resolution of the bootstrap test. See 
                 bootstrap.sample for details. (default: 0.0001). 

"""
    from bootstrap import sample
    ix = S.index.get_level_values('target').isin(inerts)
    Y = S.groupby(level='target').agg(LN_mean)
    LN_I = np.log(S.loc[ix])
    s = sample(LN_I, _all_targets_LN_mean, min_pvalue=min_pvalue)
    output = pd.Series({target:s.percentileofscore(y).drop(target, errors='ignore').mean() for target, y in Y.iteritems()}, 
                       name='LN Mean P-value')
    output.index.names = ['target']
    return 1 - output*1e-2

def percentiles(S, tiers=percentile_tiers):
    """percentiles of pandas.Series Object

    Parameter:
    ----------
    tiers : the percentiles to calculate (default: [50, 60, 70, 80, 90, 95, 99])
    """
    return pd.Series(S.quantile(q=np.array(tiers)*1e-2).values, index=pd.Index(tiers, name='Percentile').astype(str))

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

def inert_normalize(S, estimator='median', inerts=inerts):
    """Normalize tumor sizes to inert sgRNA tumor sizes.

E.g. nomralized_cells = cells.groupby(level='Mouse').transform(inert_normalize)

Variable:
---------
S : pd.Series containing Absolute Cell #s and an Index with a level entitled 
    'target'--to be used to identify tumors arising from inert sgRNAs.

Parameters:
-----------
estimator : Function to determine Central Tendency of inerts (default: 'median')

inerts : List of inert targets (default: ['Neo1', 'Neo2', 'Neo3', 'NT1']).
"""
    N = S.loc[S.index.get_level_values('target').isin(inerts)].agg(estimator)
    return S/N

