import pandas as pd
import numpy as np
from tuba_seq.tools import LN_mean, inert_normalize

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
        print("Generating", N_samples, "random samples to estimate sensitivity up to {:.4%}.".format(max_sensitivity))

    gb = ref_data.groupby(level='Mouse')
    mice = pd.Series(list(gb.groups.keys()))

    inert_sgRNAs = set(ref_data.groupby(level='target').groups.keys()) - set(active_ref_sgRNAs)
    N_active_array = np.array(N_active_sgRNAs) if hasattr(N_active_sgRNAs, '__len__') else np.array([N_active_sgRNAs])
    sorted_active = np.sort(N_active_array)
    ratios = len(active_ref_sgRNAs)/sorted_active
    
    def trial( _ ):
        df = pd.concat([gb.get_group(mouse) for mouse in mice.sample(N_mice, replace=True)])
        N_tumors = (ratios*len(df)).round().astype(int)
        max_down_sample = df.sample(N_tumors[0], replace=True)
        return pd.concat(dict(zip(sorted_active, 
            (inert_normalize(max_down_sample.iloc[0:N_t]
                                            .groupby(level='Mouse')
                                            .transform(inert_normalize, inerts=inert_sgRNAs)
                                            .groupby(level='target').agg(metric), 
                             inerts=inert_sgRNAs) 
                    for N_t in N_tumors))), names=['active_sgRNAs'])
    
    sampling_distributions = pd.DataFrame(list(map(trial, range(N_samples))))

    best_FPR_dists = sampling_distributions.T.loc[(slice(None), inert_sgRNAs), :]
    FPR_thresholds = best_FPR_dists.groupby(level='active_sgRNAs').apply(
        lambda df: df.stack().quantile(q=1-alpha/df.name/(2 if two_sided else 1)))
    return sampling_distributions.stack(level='active_sgRNAs').groupby(level='active_sgRNAs').apply(
        lambda df: (df[active_ref_sgRNAs] > FPR_thresholds[df.name]).mean())

DNA_map = {0:'A', 1:'T', 2:'C', 3:'G'}
def int_DNA_to_string(a):
    return ''.join([DNA_map[i] for i in a])

def hamming_distance(a, b):
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def longest_run(a):
    da = np.diff(a)
    iszero = np.r_[False, da == 0, False]
    dd_zero = np.diff(iszero)
    run_lengths = [1 + end - start for start, end in np.where(dd_zero)[0].reshape(-1, 2)]
    return max(run_lengths) if len(run_lengths) > 0 else 1


def create_all_sgIDs(k=8, run_length=3, min_distance=3, 
    head_end='TCCGGA', random_barcode_begin='GA', 
    restriction_enzymes={'BspEI':'TCCGGA', 'BamHI':'GGATCC'}, shuffle=False):
    from itertools import product
    import random
    """"
    Generates all possible sgIDs of length k with various characteristics. 

Generated sgIDs are: 
    1. GC-balanced (50% ATs, 50% GCs) to avoid PCR amplification bias, 
    2. Distant from all other sgIDs in the list to avoid mis-annotations from 
       read errors,
    3. Deviod of nucleotide runs (e.g. AAAA....), 
    4. Devoid of certain restriction enzyme sites.


"""
    assert k%2 == 0, "sgIDs must be an even length to have balanced GC content"
    assert min_distance >= 3, """Min distance between sgIDs must be >=3. 
This algorithm efficiently removes all sgIDs within 2 nts of each other using 
a bottom-up approach. Greater distance between sgIDs is achieved simply by
deleting any candidate sgIDs that are too proximal to ones already on the list."""
    


    all_kmers = np.array(list(product(range(4), repeat=k)))
    GC_content = (all_kmers > 1).sum(axis=1)
    balanced_kmers = all_kmers[GC_content == int(k/2)]
    run_lengths = np.array(list(map(longest_run, balanced_kmers)))
    no_run_kmers = balanced_kmers[run_lengths < run_length]
    
    if shuffle:
        random.shuffle(no_run_kmers)

    sgIDs = []
    all_neighbors = set()
    alternatives = {i:set(range(4)) - set([i]) for i in range(4)}

    contain_restriction_site = 0

    for candidate in no_run_kmers:
        sgID = int_DNA_to_string(candidate)
        extended_sgID = head_end+sgID+random_barcode_begin 
        if any((seq in extended_sgID for seq in restriction_enzymes.values())):
            contain_restriction_site += 1
            continue
        candidate_neighbors = {int_DNA_to_string(np.r_[candidate[:i], alternative, candidate[i+1:]]) for i, nuc in enumerate(candidate) for alternative in alternatives[nuc]}
        if len(all_neighbors & candidate_neighbors) == 0:
            # Candidate is not too close to existing sgID
            all_neighbors = all_neighbors | candidate_neighbors | set(sgID)
            sgIDs.append(sgID)

    drop = []

    for i, sgID in enumerate(sgIDs[:-1]):
        closest = min([hamming_distance(sgID, other) for other in sgIDs[i+1:]])
        assert closest >= 3
        if min_distance > 3:
            if closest < min_distance:
                drop.append(i)

    for i in drop[::-1]:
        sgIDs.pop(i)

    print("""For {:}mers:
{:.1%} Were balanced in GC content,
{:.1%} of the remainder had no runs of {:} or more nucleotides,
{:.1%} contained a restriction site for these enzymes: {:},
{:.1%} of the remainder avoided other sgIDs at >={:} loci. 
{:} sgIDs were identified.""".format(k, 
len(balanced_kmers)/len(all_kmers), 
len(no_run_kmers)/len(balanced_kmers), run_length,
contain_restriction_site/len(no_run_kmers), ', '.join(list(restriction_enzymes.keys())),
len(sgIDs)/(len(no_run_kmers) - contain_restriction_site), min_distance,
len(sgIDs)))

    return sgIDs

