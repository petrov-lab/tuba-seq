#!/usr/bin/env python3
import argparse, os 
import pandas as pd
import numpy as np
from tuba_seq.fastq import infer_master_read, singleMismatcher
from tuba_seq.shared import logPrint


parser = argparse.ArgumentParser(   description="""Combines DADA2 clustering outputs, annotates sgRNAs, corrects for GC-bias, performs post-clustering barcode selection, and transforms abundances to absolute cell number. All input & output data files are in CSV file format.""", 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('sgRNA_file', help='All sgRNAs used and their corresponding identifiers.')
parser.add_argument('--dir', dest='directory', default='clustered', help='Directory containing all the DADA2 clustering outputs.')
parser.add_argument('-o', '--out_file', type=str, default='combined.csv', help='CSV file with the consolidated samples and their sgID-target annotations.')
parser.add_argument('-b', '--bartender', action='store_true', help='Process Bartender (https://github.com/LaoZZZZZ/bartender-1.1) clustering output.')
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', dest='parallel', action='store_true', help='Parallelize operation.')
parser.add_argument('--indel_tolerated', type=int, default=2, help='Size of indel tolerated when attempting to annotate the sgID.')

args = parser.parse_args()
Log = logPrint(args)
csv_ext = '.csv'
if args.parallel:
    from tuba_seq.pmap import pmap as map 

if args.bartender:
    merge_rules = dict(abundance=np.sum)
    indel_tolerated = 0
else:
    merge_rules = dict(abundance=np.sum, n0=np.sum, n1=np.sum, nunq=np.sum, pval=np.prod, birth_pval=np.prod, birth_ham=np.min)
        # Dictionary of reduction operators to aggregate DADA2 clusters deemed identical
    indel_tolerated = args.indel_tolerated

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

def load_clusters_annotate_sgRNAs_and_merge(filename, columns_to_keep=list(merge_rules.keys())+['sequence']):
    """load_clusters_annotate_sgRNAs_and_merge(filename) -> dict with output & stats.

This function combines the loading, annotation, and merging steps to permit parallelization. 
"""
    df = pd.read_csv(filename, usecols=columns_to_keep) if not args.bartender else pd.read_csv(filename, usecols=[1, 3], header=0, names=['sequence', 'abundance'])
        # Immediately trim sequences down to the maximum indel tolerated. 
    master_read = infer_master_read(df['sequence'][:1000])
    start = master_read.index("N")
    stop = master_read.rindex("N")+1
    assert len(master_read[:start]) == len(master_read[stop:]), "DADA2 clusters seem to have asymmetric-lengthed flanks"
    flanking_seq_length = start
    barcode_length = stop - start
        # First, annotate all sgRNAs/barcodes that are exact matches to an sgID in the exact location
    df['ID'] = df['sequence'].str.slice(start=flanking_seq_length, stop=flanking_seq_length+sgID_length)
    df['target'] = sgID_map.loc[df['ID']].values
    df['barcode'] = df['sequence'].str.slice(start=flanking_seq_length+sgID_length, stop=flanking_seq_length+barcode_length)
        # Then, salvage the remaining tumors by permitting inexact matches and/or Indels
    failed_matches = df['target'].isnull()
    df.loc[failed_matches, :] = df.loc[failed_matches, :].apply(search_sgRNA_and_barcode, args=(flanking_seq_length,), axis=1)
    assert not df['target'].isnull().any(), "Failed to annotate an RNA"
        # Merge any tumors with the exact sgRNA annotation & random barcode sequence
    merge_rules['ID'] = lambda S: S.value_counts().argmax()
    merged = df.groupby(['target', 'barcode']).agg(merge_rules)
        # Return output & summary of merges/matching
    if args.verbose:
        print("Completed", filename.split('/')[-1], '...') 
    return dict(initial=len(df),
                clusters=merged, 
                exact_matches=len(df) - failed_matches.sum(), 
                merges=len(df) - len(merged))

Files = [f for f in os.listdir(args.directory) if csv_ext in f]
clustered_mice = map(load_clusters_annotate_sgRNAs_and_merge, [os.path.join(args.directory, f) for f in Files])

# Consolidate output based on mouse names into a single file
mice_names = [f.split(csv_ext)[0] for f in Files]
combined = pd.concat({mouse_name:output.pop('clusters') for mouse_name, output in zip(mice_names, clustered_mice)}, names=['Mouse'])
Log("Completed consolidation, sgID annotation & cluster merging.")

combined.to_csv(args.out_file+'.gz', compression='gzip')

######################### Summary Statistics ###################################
data = pd.DataFrame(clustered_mice, index=mice_names)
tallies = data.sum()
tallies['unknown_RNA'] = len(combined.loc[(slice(None), 'Unknown'), :]) 
percents = tallies/tallies.loc['initial']
Log("""{exact_matches:.2%} of DADA2 clusters perfectly matched an sgID.
{unknown_RNA:.2%} of clusters had an Unknown sgID.
{merges:.2%} of clusters were identical to another cluster after annotation & random barcode isolation. 
This number should be small (~1%.)""".format(**percents), True)

