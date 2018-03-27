#!/usr/bin/env python3
import argparse 
import pandas as pd
import numpy as np
from tuba_seq.fastq import singleMismatcher
from tuba_seq.shared import logPrint
from pathlib import Path

parser = argparse.ArgumentParser(   description="""Combines DADA2 clustering output & annotates sgRNAs.""",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('sgRNA_file', help='All sgRNAs used and their corresponding identifiers.')
parser.add_argument('--dir', dest='directory', type=Path, default='clustered', help='Directory containing all the DADA2 clustering outputs.')
parser.add_argument('-o', '--out_file', type=str, default='combined.csv', help='CSV file with the consolidated samples and their sgID-target annotations.')
parser.add_argument('-i', '--input', default='dada2', choices=['dada2', 'bartender', 'derep'], help='Barcode output to handle.')
parser.add_argument('-d', '--derep', action='store_true', help="Process derep'ed .rds files that have gone through preprocessing.")
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', dest='parallel', action='store_true', help='Parallelize operation.')
parser.add_argument('--indel_tolerated', type=int, default=2, help='Size of indel tolerated when attempting to annotate the sgID.')
#parser.add_argument('--master_read', type=str, default='', help='Read outline to expect.')
parser.add_argument('--flank', type=int, default=4, help='Expected beginning and end length of reads.')

args = parser.parse_args()
Log = logPrint(args)
if args.parallel:
    from tuba_seq.pmap import pmap as map 

if args.input == 'bartender':
    merge_rules = dict(abundance=np.sum)
    indel_tolerated = 0
    def read_input(filename): 
        return pd.read_csv(filename, usecols=[1, 3], header=0, names=['sequence', 'abundance'])
    file_glob = '*.csv*'
elif args.input == 'dada2':
    merge_rules = dict(abundance=np.sum, n0=np.sum, n1=np.sum, nunq=np.sum, pval=np.prod, birth_pval=np.prod, birth_ham=np.min)
        # Dictionary of reduction operators to aggregate DADA2 clusters deemed identical
    indel_tolerated = args.indel_tolerated
    def read_input(filename):
        return pd.read_csv(filename, usecols=list(merge_rules.keys()+['sequence']))
    file_glob = '*.csv*'
elif args.input == 'derep':
    merge_rules = dict(abundance=np.sum)
    indel_tolerated = args.indel_tolerated
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    readRDS = robjects.r['readRDS']
    names = robjects.r['names']
    def read_input(filename):
        robj = readRDS(str(filename))[0]
        return pd.DataFrame(dict(   sequence=pandas2ri.ri2py(names(robj)),
                                    abundance=pandas2ri.ri2py(robj)))
    file_glob = '*.rds'

sg_info = pd.read_csv(args.sgRNA_file)
sg_info['ID'] = sg_info['ID'].str.upper()

sgID_lengths = sg_info['ID'].str.len()
sgID_length = sgID_lengths.values[0]
assert (sgID_lengths == sgID_length).all(), 'This program expects all sgIDs to be the same length.'

sgID_map = sg_info.set_index('ID')['target']

sg_matchers = pd.Series({target:singleMismatcher(sg_id) for sg_id, target in sgID_map.items()})

def search_sgRNA_and_barcode(row, flanking_seq_length):
    """search_sgRNA_and_barcode(DADA2_cluster_row) -> row with correct sgRNA/barcode

Identifies the sgRNA corresponding to the various sgIDs enumerated in the sgRNA
file. Tolerates a single mismatch and small indel when searching for the sgID. 
"""
    DNA = row["sequence"][flanking_seq_length-indel_tolerated:]
    positions = sg_matchers.apply(lambda func: func.find(DNA[:sgID_length+2*indel_tolerated+1]))
    finds = positions[positions >= 0]
    exact_locations = finds[finds == indel_tolerated]
    if len(exact_locations) == 1:   # sgID is a single mutation in the expected location; barcode unchanged
        row['target'] = exact_locations.index[0]
    elif len(finds) == 1:           # sgID is not in the expected location, but there is only 1 reasonable match
                                    # These two cases are separated because you could have 2 reasonable matches, but only 1 in the expected location
        row['target'] = finds.index[0]
        start = finds[0]
        row['ID'] = DNA[start:start+sgID_length]
        barcode_length = len(DNA) - 2*flanking_seq_length
        row['barcode'] = DNA[start+sgID_length:start+barcode_length]
    else:                           # There are 0 or 2+ reasonable, yet imperfect sgID matches.
        row['target'] = "Unknown" 
    return row

def infer_master_read(DNAs, max_random_base_frequency=0.667):
    counts = pd.DataFrame({i:DNAs.str.get(i).value_counts() for i in range(len(DNAs[0]))})
    PWM = counts/counts.loc[['A', 'C', 'G', 'T']].sum()
    master_read = ''.join(PWM.apply(lambda col: col.idxmax() if col.max() > max_random_base_frequency else 'N'))
    for ID in sg_info['ID']:
        if ID in master_read: 
            begin, ID, end = master_read.partition(ID)
            master_read = begin + len(ID)*'N' + end
    assert len(master_read.partition('N')[0]) == len(master_read.rpartition('N')[2]), "Inferred master read has uneven flanks:\n"+master_read+PWM.to_string()
    return master_read

def load_clusters_annotate_sgRNAs_and_merge(filename):
    """load_clusters_annotate_sgRNAs_and_merge(filename) -> dict with output & stats.

This function combines the loading, annotation, and merging steps to permit parallelization. 
"""
    sample = filename.stem
    df = read_input(filename)
        # Immediately trim sequences down to the maximum indel tolerated. 
    #master_read = infer_master_read(df['sequence'][:10000])
    start = args.flank#master_read.index("N")
    stop = len(df['sequence'].iloc[0]) - args.flank #master_read.rindex("N")+1
    flanking_seq_length = start
    barcode_length = stop - start
        # First, annotate all sgRNAs/barcodes that are exact matches to an sgID in the exact location
    df['ID'] = df['sequence'].str.slice(start=flanking_seq_length, stop=flanking_seq_length+sgID_length)
    #df['target'] = sgID_map.loc[df['ID']].values
    df['target'] = sgID_map.reindex(df['ID']).values
    df['barcode'] = df['sequence'].str.slice(start=flanking_seq_length+sgID_length, stop=flanking_seq_length+barcode_length)
        # Then, salvage the remaining tumors by permitting inexact matches and/or Indels
    failed_matches = df['target'].isnull()
    df.loc[failed_matches, :] = df.loc[failed_matches, :].apply(search_sgRNA_and_barcode, args=(flanking_seq_length,), axis=1)
    assert not df['target'].isnull().any(), "Failed to annotate an RNA"
        # Merge any tumors with the exact sgRNA annotation & random barcode sequence
    merge_rules['ID'] = lambda S: S.value_counts().idxmax()
    merged = df.groupby(['target', 'barcode']).agg(merge_rules)
        # Return output & summary of merges/matching
    if args.verbose:
        print("Completed", sample, '...') 
    return dict(Sample=sample,
                initial=len(df),
                clusters=merged, 
                exact_matches=len(df) - failed_matches.sum(), 
                merges=len(df) - len(merged)) #,
#                master_read=master_read)

Files = list(args.directory.glob(file_glob))
clustered_samples = map(load_clusters_annotate_sgRNAs_and_merge, Files)

# Consolidate output based on sample names into a single file
combined = pd.concat({dic['Sample']:dic.pop('clusters') for dic in clustered_samples}, names=['Sample'])
Log("Completed consolidation, sgID annotation & cluster merging.")

combined.to_csv(args.out_file+'.gz', compression='gzip')

# Check Master Read Inferences
data = pd.DataFrame(clustered_samples).set_index("Sample")
#master_reads = data.pop('master_read').str.replace("N", '.')
#MR_counts = master_reads.value_counts()
#if len(MR_counts) == 1:
#    Log("All master reads matched ({:})".format(MR_counts.values[0]))
#else:
#    from warnings import warn
#    main_read = MR_counts.idxmax()
#    s = "Master reads did not all match! Used {:} {:} times, and observed the following exceptions:\n".format(main_read, MR_counts.max())
#    s += master_reads.loc[master_reads != main_read].to_string()
#    warn(s)
#    Log(s)


######################### Summary Statistics ###################################
tallies = data.sum()
tallies['unknown_RNA'] = len(combined.loc[(slice(None), 'Unknown'), :]) 
percents = tallies/tallies.loc['initial']
Log("""{exact_matches:.2%} of DADA2 clusters perfectly matched an sgID.
{unknown_RNA:.2%} of clusters had an Unknown sgID.
{merges:.2%} of clusters were identical to another cluster after annotation & random barcode isolation. 
This number should be small (~1%.)""".format(**percents), True)

