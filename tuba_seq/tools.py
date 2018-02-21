import numpy as np
import pandas as pd
import io, sys, collections
from tuba_seq.functions import LN_Mean_P_values, percentiles, LN_mean, inerts 

def LN_Mean_Summary(S, inerts=inerts, min_FWER=0.0001, correction='Bonferroni'):
    """Hypothesis test of increased growth for sgRNAs within the screen, assuming 
a Lognormal distribution of tumor sizes. 

This function is similar to tuba_seq.functions.LN_Mean_P_values, however it 
corrects for multiple-hypothesis testing, and generates a more descriptive table. 
This function uses the bootstrap sampling technique to estimate a null sampling 
distribution from the sgInert tumors sizes. 

Note: This hypothesis test was not used in Rogers et al. 2017. This test 
*assumes* that sgInert tumors represent a null sampling distribution and, 
therefore, cannot test whether any sgInert guide deviates in growth from the 
other inerts. 
    
This test tends to be more conservative in P-value estimates because it attempts
to account for the typical amount of off-target endonuclease activity insofar as 
it alters cell fitness. For example, sgNeo2 tumors tend to grow slightly slower 
than other inerts. By estimating the mean size of sgNeo2 tumors separate from 
the other inerts (and by generating a sampling distribution from all sgInerts 
separately), this unintended change in fitness is incorporated into the null 
sampling distribution, widening the null distribution and, thus, reducing the 
likelihood that target sgRNAs deviate significantly from the null distribution. 

This test (not the Nat Methods test) is recommended. However, if one of the 
`inerts` within the screen changes tumor growth inordinately, then you will not 
observe any statistically-significant growth effects. So if you are using a new 
sgInert, you may want to first excluded it from `inerts` to affirm that it 
minimally alters tumor growth. 
    
Parameters:
-----------

inerts : sgRNAs that are not expected to alter tumor sizes to be used as the 
    null sampling distribution (default: ['Neo1-3', 'NT1', 'NT3'])

min_FWER : Determines number of bootstrap samples to draw. Family-Wise Error Rate
    (FWER) estimates below min_FWER are not reliably--more bootstrapping samples 
    need to be drawn. The default resolution can take hours to run even on a 
    modern CPU architecture. See bootstrap.sample for details. (default: 0.0001). 

correction : Correction method for FWER (default: 'Bonferroni', other options: 
    'Sidak', 'Holm', 'Hochberg'. See tuba_seq.bootstrap for details. The # of
    hypotheses is estimated from the non-inert sgRNAs in the input array.
"""
   out = S.groupby(level='target').agg(LN_mean)
   from 
   min_pval = #####
   df = pd.concat({
        'LN Mean (absolute cell no.)':out, 
        'LN Mean (Relative to sgInerts)':out.groupby(level='target').transform(inert_normalize)},
        'One-Sided raw P-value': LN_mean_P_values(S, inert=inerts, min_pvalue=min_pvalue))



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


