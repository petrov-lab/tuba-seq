from numpy cimport ndarray 
from tuba_seq.pmap import pmap as map 
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

def percentiles(S, tiers=percentile_tiers):
    """percentiles of pandas.Series Object

    Parameter:
    ----------
    tiers : the percentiles to calculate (default: [50, 60, 70, 80, 90, 95, 99])
    """
    return pd.Series(S.quantile(q=np.array(tiers)*1e-2).values, index=pd.Index(tiers, name='Percentile').astype(str))

def best_power_law_fit(S, resolution=400, percentile_min=10, percentile_max=99, sigma_threshold=0.05):
    """Returns best-fitting powerlaw.Distribution object of the data.

    The built-in optimization routine for the powerlaw package often fails to find
    the best/most-reasonable fit. This function finds an optimal fit via brute-force
    and restricts fits to sane parameters for TuBa-seq purposes. 

    Parameter:
    ----------
    resolution : The number of evenly-spaced values for `xmin` to test. (default: 400)

    percentile_min : The smallest `xmin` value to consider, as a percentile of the
        distribution (default: 10)

    percentile_max: The largest `xmin` percentile to consider (default: 99)

    sigma_threshold : Upper limit on the standard error of the power law fit. 
        (default: 0.05)
    """
    import powerlaw as pl
    X_lims = percentiles(S, tiers=[percentile_min, percentile_max])
    grid = np.linspace(*X_lims, num=resolution)
    Fits = [pl.Fit(S, xmin=xmin, sigma_threshold=sigma_threshold, fit_method='Likelihood') for xmin in grid]
    log_likelihoods = np.array([fit.distribution_compare('power_law', 'lognormal', normalized_ratio=True)[0] for fit in Fits])
    return Fits[log_likelihoods.argmax()], log_likelihoods.max()

def best_parametric_distribution(S, consider=['norm', 'expon', 'logistic', 'gumbel', 'lognorm', 'loglogistic']):
    """Anderson-Darling test to identify the best parametric distribution of the data.

    Parameter:
    ----------
    consider : List of distributions to compare. Will return best fit.
    """
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

def inert_normalize(S, estimator='mean', inerts=inerts):
    N = S.loc[S.index.get_level_values('target').isin(inerts)].agg(estimator)
    return S/N

def load_tumors(filename='tumor_sizes.csv.gz', metadata_filename='sample_metadata.csv', 
                meta_columns=[], drop=None, merger=None, merge_func='mean'):

    tumors = pd.read_csv(filename)
    meta_df = pd.read_csv(metadata_filename).set_index("Sample").sort_index()
   
    ix_names = ['Sample', 'target', 'barcode']
    if drop is not None:
        if type(drop) == str:
            if not drop in meta_df.columns:
                raise LookupError("drop `{:}` is not a column in the metadata file.".format(drop))
            drop = meta_df.query('@drop or @drop == "y"').index
        tumors = tumors.loc[tumors.isin(drop)]

    if merger is not None:
        if type(merger) == str:
            if not merger in meta_df.columns:
                raise LookupError("merger `{:}` is not a column in the metadata file.".format(merger))
            merger = meta_df[merger]
        elif merger.name is None:
            raise ValueError("Merger must have a `name` attribute to label this new column.")
        
        tumors.insert(0, merger.name, merger.loc[tumors['Sample']].values)
        meta_df = meta_df.reset_index().set_index(merger.name)
        ix_names[0] = merger.name
        tumors = tumors.groupby(level=ix_names).agg(merge_func)

    if meta_columns != []:
        tumors = pd.concat([meta_df.loc[tumors[ix_names[0]], meta_columns], tumors], axis=1)
        ix_names = meta_columns + ix_names
    return tumors.set_index(ix_names)['Cells'].sort_index()


PERMISSIBLE_UNCERTAINTY = 0.21

def power_analysis(ref_data, active_ref_sgRNAs, N_mice, N_active_sgRNAs, 
                   metric=LN_mean, alpha=0.05, two_sided=True, max_sensitivity=0.99, N_samples=None):
    """Sensitivity of Tuba-seq to proposed experiment. 
    
    Returns the Sensitivity (TPR) of a hypothetical Tuba-seq experiment by down-
    sampling tumor sizes from a reference dataset. The effect sizes are the true 
    effect sizes of the active sgRNAs within the reference dataset. If you want to
    project results for a hypothetical effect size, I recommend multiplying the active
    sgRNAs by constants that impart the effect size of interest. To reduce the number
    of permutation tests, this analysis does not model the bootstrapping method that
    is used to naively identify p-values in a real Tuba-seq dataset. Instead, it 
    models many hypothetical experiments to generate a sampling distribution of inert
    sgRNAs that defines the size threshold for the desired FPR, and also a sampling 
    distribution of each active sgRNA to identify the sensitivity above this FPR. As
    such, it may be a little more conservative, but also should more accurately capture
    the statistical noise introduced via off-target sgRNA cutting insofar as it is
    captured by variation in inert sgRNA mean sizes. 

    Variables:
    ----------
    ref_data : pandas.Series of tumor sizes indexed by `Mouse` and `target` genes.
    
    active_ref_sgRNAs : Iterable of active sgRNAs in the reference dataset.

    N_mice : Number of mice in hypothetical experiment.

    N_active_sgRNAs : Number of active sgRNAs in hypothetical experiment. 

    Parameters:
    -----------
    metric : Summary statistical measure to define growth (default: LN_mean)

    alpha : Desired Specificity (FPR) of survey (default: 0.05)

    two_sided : Use two-sided statistical test for increased growth (default: True)

    N_samples : Number of random samples to generate for determining Sensitivity. 
        (default: None -- See below).

    max_sensitivity : This parameter is only relevant when N_samples == None. 
        Uses a common heuristic to define the maximum trust-able Sensitivity level.
        Projecting higher TPRs requires more down-samplings to adequately model the 
        left-most tail of the sampling distribution (False Negatives) from which the 
        Sensitivity rate is determined. (default: 0.99, I.e. you can trust Sensitivity
        values reported up to 99% sensitivity, but if you wanted to know where your
        Tuba-seq experiment becomes 99.9% sensitive, then you must generate more
        random samplings.)
    """

    if N_samples is None:
        N_samples = int(np.ceil(4/((1-max_sensitivity)*PERMISSIBLE_UNCERTAINTY**2)))
        print("Running", N_samples, "Samples")

    gb = ref_data.groupby(level='Mouse')
    mice = pd.Series(list(gb.groups.keys()))

    inerts = set(ref_data.groupby(level='target').groups.keys()) - set(active_ref_sgRNAs)
    
    ratio = len(active_ref_sgRNAs)/N_active_sgRNAs
    def analyze(df):
        down_sample = df.sample(int(round(ratio*len(df))), replace=True)
        return inert_normalize(down_sample.groupby(level='target').agg(metric), inerts=inerts)

    Y = pd.DataFrame([analyze(pd.concat([gb.get_group(mouse) for mouse in mice.sample(N_mice, replace=True)])) for _ in range(N_samples)])
    
    alpha_B = (alpha/N_active_sgRNAs/2) if two_sided else (alpha/N_active_sgRNAs)
    FPR_threshold = Y[list(inerts)].stack().quantile(q=1-alpha_B) 
    print('\n', N_active_sgRNAs, FPR_threshold, alpha_B, Y[list(inerts)].stack().quantile(q=0.99))
    return (Y[active_ref_sgRNAs] > FPR_threshold).mean()
