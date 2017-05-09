#!/usr/bin/python3
import argparse, os
import pandas as pd
import numpy as np
from tuba_seq.fastq import singleMismatcher
from params import barcode_length, merge_rules 

tuba_seq_dir = os.path.dirname(__file__)

parser = argparse.ArgumentParser(   description="""Combines DADA2 clustering outputs and annotates sgRNAs. All input & output data files should be in CSV file format.""", 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir', dest='directory', default='clustered', help='Directory containing all the DADA2 clustering outputs.')
parser.add_argument('--indel_tolerated', type=int, default=2, help='Size of indel tolerated when attempting to annotate the sgID.')
parser.add_argument('-o', '--out_file', type=str, default='combined.csv')
parser.add_argument('--sgRNA_file', default='{:}/sgRNA_info.csv'.format(tuba_seq_dir), help='All sgRNAs used and their corresponding identifiers.')
parser.add_argument('-p', '--parallel', dest='parallel', action='store_true', help='Parallelize operation.')
parser.add_argument('-b', '--bartender', action='store_true', help='Process Bartender (https://github.com/LaoZZZZZ/bartender-1.1) clustering output.')
args = parser.parse_args()

csv_extension = '.csv'

if args.parallel:
    from tuba_seq.pmap import pmap as map 

sg_info = pd.read_csv(args.sgRNA_file)
sg_info['ID'] = sg_info['ID'].str.upper()

sgID_lengths = sg_info['ID'].str.len()
sgID_length = sgID_lengths.values[0]
assert (sgID_lengths == sgID_length).all(), 'This program expects all sgIDs to be the same length.'

sg_info['RNA'] = sg_info['sgRNA'].str.replace('sg', '')
sgRNA_map = sg_info.set_index('ID')['RNA']

sg_matchers = [singleMismatcher(sg_id) for sg_id in sgRNA_map.index]

correct_sgIDs = frozenset(sg_info['ID'])
indel_tolerated = args.indel_tolerated if not args.bartender else 0

def search_sgRNA_and_barcode(row):
    """identify_sgRNA_and_barcode(DADA2_cluster_row) -> row with correct sgRNA/barcode

Identifies the sgRNA corresponding to the various sgIDs enumerated in the sgRNA
file. Tolerates a single mismatch and small indel when searching for the sgID. 
"""
    DNA = row["sequence"]
    positions = np.array([sg_matcher.find(DNA[0:sgID_length+2*indel_tolerated+1]) for sg_matcher in sg_matchers])
    exact_location_sgRNAs = sgRNA_map.ix[positions == indel_tolerated]
    if len(exact_location_sgRNAs) == 1: # sgID is a single mutation in the expected location; barcode unchanged
        row['RNA'] = exact_location_sgRNAs.values[0]
        return row
    finds = np.where(positions >= 0)[0]
    if len(finds) == 1:             # sgID is not in the expected location, but there is only 1 reasonable match
                                    # These two cases are separated because you could have 2 reasonable matches, but only 1 in the expected location
        position_ix = finds[0]
        position = positions[position_ix]
        row['RNA'] = sgRNA_map.iloc[position_ix]
        row['barcode'] = DNA[position+sgID_length:position+barcode_length]
    else:                           # There are 0 or 2+ reasonable, yet imperfect sgID matches.
        row['RNA'] = "Unknown" 
    return row

def load_clusters_annotate_sgRNAs_and_merge(filename, usecols=list(merge_rules.keys())+['sequence']):
    """load_clusters_annotate_sgRNAs_and_merge(filename) -> dict with output & stats.

This function combines the loading, annotation, and merging steps to permit parallelization. 
"""
    df = pd.read_csv(filename, usecols=usecols if not args.bartender else None)
    if args.bartender:
        df.columns = ['Cluster.ID', 'sequence', 'Cluster.score', 'abundance']
        merge_rules = dict(abundance=np.sum)
        df.pop("Cluster.score")
    num_clusters = len(df)
    seq_length = len(df['sequence'].values[0]) 
    flanking_seq_length = int( (seq_length - barcode_length)/2 ) if not args.bartender else 0
    df['sequence'] = df['sequence'].str.slice(start=flanking_seq_length - indel_tolerated, stop=flanking_seq_length + barcode_length + indel_tolerated)

    df['RNA'] = sgRNA_map.loc[df['sequence'].str.slice(start=indel_tolerated, stop=indel_tolerated+sgID_length)].values
    df['barcode'] = df['sequence'].str.slice(start=indel_tolerated+sgID_length, stop=indel_tolerated+barcode_length)
    failed_matches = df['RNA'].isnull()
    num_failed_matches = failed_matches.sum()
    df.loc[failed_matches, :] = df.loc[failed_matches, :].apply(search_sgRNA_and_barcode, axis=1)
    assert not df['RNA'].isnull().any(), "Failed to annotate an RNA"
    
    df.set_index(['RNA', 'barcode'], inplace=True)
    merged = df.groupby(level=['RNA', 'barcode']).agg(merge_rules)
    if not args.bartender:
        merged['n10'] = merged['n0'] + merged.pop('n1')
    
    return dict(clusters=merged, 
                exact_matches=num_clusters - num_failed_matches, 
                merged_clusters=num_clusters - len(merged))

Files = [f for f in os.listdir(args.directory) if csv_extension in f]
clustered_mice = map(load_clusters_annotate_sgRNAs_and_merge, [args.directory+'/'+f for f in Files])

mice_names = [f.split(csv_extension)[0].split('_')[0] for f in Files]

combined = pd.concat({mouse_name:output['clusters'] for mouse_name, output in zip(mice_names, clustered_mice)}, names=['Mouse'])

all_merges = sum((output['merged_clusters'] for output in clustered_mice))
all_exact_matches = sum((output['exact_matches'] for output in clustered_mice))
unknown_RNAs = len(combined.loc[(slice(None), 'Unknown'), :])
total_clusters = len(combined) + all_merges

print("""{:.2%} of DADA2 clusters perfectly matched an sgID.
{:.2%} of clusters had an Unknown sgID.
{:.2%} of clusters were identical to another cluster after their sgRNA & barcode were annotated & isolated. This number should be small (i.e. ~1%.)""".format(
all_exact_matches/total_clusters,
unknown_RNAs/total_clusters,
all_merges/total_clusters))

combined.to_csv(args.out_file+'.gz', compression='gzip')

