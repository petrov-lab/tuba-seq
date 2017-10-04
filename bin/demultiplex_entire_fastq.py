#!/usr/bin/env python3

import os, argparse
from tuba_seq.fastq import singleMismatcher 
from tuba_seq.shared import logPrint
import pandas as pd

parser = argparse.ArgumentParser(description="Split paired-end read files by Illumina indecies.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("forward_read_file", help='FASTQ file of forward reads')
parser.add_argument("reverse_read_file", help='FASTQ file of reverse reads')
parser.add_argument('barcode_file', type=str, help='Tab-delimited file with sample_name, barcode pairs.')
parser.add_argument("forward_read_dir", default='forward_reads', help='Directory to put split forward reads.')
parser.add_argument("reverse_read_dir", default='reverse_reads', help='Directory to put split reverse reads.')
###############################################################################

args = parser.parse_args()

Log = logPrint(args)

samples = pd.read_csv(args.barcode_file, sep='\t', names=['Samples', 'Index'], index_col=1)['Samples']
length = max(samples.index.str.len())
sample_matchers = {sample:singleMismatcher(index) for sample, index in samples.items()}

Log("Inferred length of barcodes to be {:} nts.".format(length))
forward_ext = args.forward_reads.split('.')[-1]
reverse_ext = args.reverse_reads.split('.')[-1]

assert forward_ext == reverse_ext, 'compression/filetype of forward reads must match reverse reads'

fastq_ext = '.fastq'
if forward_ext =='gz':
        from gzip import open       # We don't need to specify binary mode...it will automatically assume.    
        fastq_ext += '.gz'
elif forward_ext == 'bz2':
        from bz2 import open
        fastq_ext += '.bz2'

counts = dict(mismatches=0, errors=0)
def identify_barcode(bc1, bc2):
    if bc1 == bc2: 
        bc = bc1
    else:
        if 'N' in bc1:
            bc1 = ''.join([r_i if f_i == "N" else f_i for f_i, r_i in zip(bc1, bc2)])
        if 'N' in bc2:
            bc2 = ''.join([f_i if r_i == "N" else r_i for f_i, r_i in zip(bc1, bc2)])
        if bc1 == bc2:  # N-resurrectable
            bc = bc1
        else:
            counts['mismatches'] += 1
            return ('conflict' if bc2 in samples else samples[bc1]) if bc1 in samples else samples.get(bc2, 'unknown') 
    if bc in samples:
        return samples[bc]
    else:
        matches = {sample for sample, matcher in sample_matchers.items() if matcher.find(bc) != -1}
        if len(matches) == 1:
            counts['errors'] +=1 
            return matches.pop()
        else:
            return 'unknown'

extras = ['conflict', 'unknown']
filenames = samples.tolist()+extras

file_tallies = {filename:0 for filename in filenames}

os.makedirs(args.forward_read_dir, exist_ok=True)
os.makedirs(args.reverse_read_dir, exist_ok=True)
forward_files = {filename:open(os.path.join(args.forward_read_dir, filename)+fastq_ext, 'wb') for filename in filenames}
reverse_files = {filename:open(os.path.join(args.reverse_read_dir, filename)+fastq_ext, 'wb') for filename in filenames}
with open(args.forward_reads) as forward_file, open(args.reverse_reads) as reverse_file:
    for i, (forward_line, reverse_line) in enumerate(zip(forward_file, reverse_file)):
        if i%4 == 0:
            filename = identify_barcode(forward_line.decode('ascii').split(':')[-1][:length], reverse_line.decode('ascii').split(':')[-1][:length])
            file_tallies[filename] += 1
        forward_files[filename].write(forward_line)
        reverse_files[filename].write(reverse_line)

for f in forward_files.values():
    f.close()
for f in reverse_files.values():
    f.close()

S = pd.Series(file_tallies)

total = S.sum()
successes = total - S[extras].sum()

Log("""Split {:,} reads ({:.2%} successfully).')
{:.2%} of forward-reverse barcodes matched.")
"{:.2%} of barcodes were saved by tolerating single-nucleotide errors in the barcode. 
Reads were split into the following files:""".format(
total, successes/total, 1 - counts['mismatches']/total, counts['errors']/total ))
percentages = S/total
for k, v in percentages.items():
    Log('{:<20} {:.2%}'.format(k+':', v))


