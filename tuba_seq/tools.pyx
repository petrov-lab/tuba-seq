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

def xformed_LN_mean(ndarray[dtype=double, ndim=1] LN_x):
    """deprecated"""
    cdef:
        double L = len(LN_x)
        ndarray[dtype=double, ndim=1] One = np.ones_like(LN_x)
    X = LN_x.dot(One)/L
    X2 = LN_x.dot(LN_x)/L - X*X
    return np.exp(X + 0.5*X2)

def LN_mean_P_values(S, inerts=inerts, min_pvalue=0.0001):
    from bootstrap import sample
    ix = S.index.get_level_values('target').isin(inerts)
    Y = S.groupby(level='target').agg(LN_mean)
    LN_I = np.log(S.loc[ix])
    def all_targets_LN_mean(ln_S):
        return ln_S.groupby(level='target').agg(lambda S: xformed_LN_mean(S.values))
    s = sample(LN_I, all_targets_LN_mean, min_pvalue=min_pvalue)
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

def best_parametric_distribution(S, consider=['norm', 'expon', 'logistic', 'gumbel', 'extreme1', 'lognorm', 'loglogistic']):
    """Anderson-Darling test to identify the best parametric distribution of the data.

    Parameter:
    ----------
    consider : List of distributions to compare. Will return best fit.
    """
    from scipy.stats import anderson
    log_distributions = {'lognorm', 'loglogistic'}
    fits = {dist:anderson(np.array(np.log(S) if dist in log_distributions else S), dist=dist[(3 if dist in log_distributions else 0):]).statistic for dist in consider}
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

def load_tumors(filename='tumor_sizes.csv.gz', metadata='sample_metadata.csv', 
                meta_columns=[], drop=None, merger=None, merge_func='mean'):
    """Loads & annotates tumor size datafile, returns pd.Series of Absolute Cell #.

Parameters:
-----------
filename : CSV file containing all tumor sizes and their sample/sgRNA target/
    barcodes (default: "tumor_sizes.csv.gz'--default of final_processing.py).

metadata : CSV file containing sample metadata (default: 'sample_metadata.csv').

meta_columns : Columns in the metadata file to append to the Multi-Index of the 
    output pd.Series. For example, I typically include a 'Genotype' and 'Time' 
    (of tumor growth in weeks) column in the metadata file, which I use to group
    and slice tumor sizes in various analyses (default: []).

drop : Samples to exclude from output. Either column name in metadata file (with 
    boolean or y/n entries, or container of Sample names. (default: None).

merger : Mergers replicate samples. Either column name in metadata file (with
    final sample names as entries, or dictionary with Sample -> merged_name
    mappings. (default: None).

merge_func : Function to merge tumors in replicate samples (default: 'mean').
"""
    tumors = pd.read_csv(filename)
    meta_df = pd.read_csv(metadata).set_index("Sample").sort_index()
    ix_names = ['Sample', 'target', 'barcode']
    if drop is not None:
        if type(drop) == str:
            if not drop in meta_df.columns:
                raise LookupError("drop `{:}` is not a column in the metadata file.".format(drop))
            drop = meta_df.query('@drop or @drop == "y"').index
        tumors = tumors.loc[tumors.isin(drop)]
    
    if meta_columns != []:
        tumors = tumors.join(meta_df.loc[tumors['Sample'], meta_columns].reset_index(drop=True))
    
    if merger is not None:
        if type(merger) == str:
            if not merger in meta_df.columns:
                raise LookupError("merger `{:}` is not a column in the metadata file.".format(merger))
            merger = meta_df[merger]
        elif merger.name is None:
            raise ValueError("Merger must have a `name` attribute to label this new column.")
        
        tumors.insert(0, merger.name, merger.loc[tumors['Sample']].values)
        return tumors.groupby(meta_columns+[merger.name]+ix_names[1:])["Cells"].agg(merge_func)
    else:     
        return tumors.set_index(meta_columns+ix_names)['Cells'].sort_index()


