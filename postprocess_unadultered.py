#!/usr/bin/python3
import argparse, os
import pandas as pd
import numpy as np
from fastq import singleMismatcher
from params import barcode, sgID

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Combines DADA2 clustering outputs and annotates sgRNAs.')
parser.add_argument('--ext', dest='extension', argument_default='.csv')
parser.add_argument('--dir', dest='directory', argument_default='clustered')
parser.add_argument('--flank', dest='cluster_flank', type=int, argument_default=7)
parser.add_argument('--out', dest='out_file', argument_default='combined.csv')
parser.add_argument('--sgRNAs', dest='sgRNA_file', argument_default='sgRNA_info.csv', help='CSV file containing ')
parser.add_argument('--parallel', dest='parallel', action='store_true', help='Parallelize operation (requires multiprocessing module.)')

args = parser.parse_args()

if args.parallel


merge_rules = dict(abundance=np.sum, n0=np.sum, n1=np.sum, nunq=np.sum, pval=np.prod, birth_pval=np.prod, birth_ham=np.min)

sg_info = pd.read_csv(args.sgRNA_file)
sg_info['ID'] = sg_info['ID'].str.upper()
sg_info['RNA'] = sg_info['sgRNA'].str.replace('sg', '')
sgRNA_map = sg_info.set_index('ID')['RNA']

sg_matchers = [singleMismatcher(sg_id) for sg_id in sgRNA_map.index]

sgID_length = len(params.sgID)
barcode_length = len(params.barcode)

correct_sgIDs = frozenset(sg_info['ID'])

def identify_RNA_barcode(DNA, cluster_flank, indel_tolerated=2):
    putative_sgID = DNA[cluster_flank:cluster_flank+sgID_length]
    putative_barcode = DNA[cluster_flank+sgID_length:cluster_flank+barcode_length]
    
    if putative_sgID in correctIDs:       #sgID is a legit sequence in the correct place
        return sgRNA_map[guess_1], guess_barcode 
    guess_2 = DNA[cluster_flank - indel_tolerated: cluster_flank + L_sgID + indel_tolerated]
    positions = np.array([sg_matcher.find(guess_2) for sg_matcher in sg_matchers])
    all_matches = positions[positions >= 0]
    top_matches = all_matches[positions[all_matches] == indel_tolerated]
    if len(top_matches) == 1:       #sgID is a single mutation in the correct place
        return sgRNA_map.ix[top_matches[0]], guess_barcode 
    elif len(all_matches == 1):     #sgID is not in the correct place, but there is exactly 1 reasonable sgID within a 2 nt vicinity.
        match = all_matches[0]
        position = positions[match] + cluster_flank - indel
        return sgRNA_map.ix[match], DNA[position+L_sgID:position+degen_len]
    return "Unknown", guess_barcode #There is not a single reasonable sgIDs in this sequence (maybe 0, maybe >=2). 











usecols = list(merge_rules.keys())

def load_and_annotate_sgRNAs(filename):
    usecols.append('sequence')
    df = pd.read_csv(filename, index_col=0, usecols=usecols)
    RNAs_and_barcodes = df['sequence'].apply(identify_RNA_barcode, args=(args.cluster_flank,))
    df.index = pd.MultiIndex.from_tuples(RNAs_and_barcodes, names=['RNA', 'barcode'])
    return df

Files = [f for f in os.listdir(args.directory) if args.extension in f]
clustered_mice = map(load_and_annotate_sgRNAs, [args.directory+'/'+f for f in Files])

mice_names = [f.split(args.extension)[0].split('_')[0] for f in Files]

combined = pd.concat({mouse_name:df for mouse_name, df in zip(mice_names, clustered_mice)}, names=['Mouse'])

merged = combined.groupby(level=['Mouse', 'RNA', 'barcode']).agg(merge_rules)
print("{:.2%} of DADA2 clusters were identical after the sgRNA was annotated and the random nucleotides were isolated. This number should be <1%.".format(1 - len(merged)/len(combined)))

merged['n10'] = merged['n0'] + merged.pop('n1')

combined.to_csv(args.out_file+'.gz', compression='gzip')

