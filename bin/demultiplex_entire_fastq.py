#!/usr/bin/env python3

import os, argparse
from tuba_seq.fastq import singleMismatcher 
from tuba_seq.shared import smart_open, logPrint
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser(description="Split paired-end read files by Illumina indecies.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("forward_read_file", help='FASTQ file of forward reads')
parser.add_argument("reverse_read_file", help='FASTQ file of reverse reads')
parser.add_argument('barcode_file', type=str, help='Tab-delimited file with sample_name, barcode pairs.')
parser.add_argument("--forward_read_dir", default='forward_reads', help='Directory to put split forward reads.')
parser.add_argument("--reverse_read_dir", default='reverse_reads', help='Directory to put split reverse reads.')
parser.add_argument('--compression', default='gz', choices=['bz2', 'gz', 'lzma', 'none'], help='Compression algorithm for output.')
###############################################################################

args = parser.parse_args()
Log = logPrint(args)
samples = pd.read_csv(args.barcode_file, sep='\t', names=['Samples', 'Index'], index_col=1)['Samples']
length = max(samples.index.str.len())
sample_matchers = {sample:singleMismatcher(index) for sample, index in samples.items()}

Log("Inferred length of barcodes to be {:} nts.".format(length))

mismatches, errors = 0, 0

def identify_barcode(bc1, bc2):
    global mismatches, errors

    if bc1 == bc2: 
        bc = bc1
    else:
            # Resolve unknown "N" bases
        if 'N' in bc1:
            bc1 = ''.join([r_i if f_i == "N" else f_i for f_i, r_i in zip(bc1, bc2)])
        if 'N' in bc2:
            bc2 = ''.join([f_i if r_i == "N" else r_i for f_i, r_i in zip(bc1, bc2)])
        if bc1 == bc2:  
            # Conflict resolved by ignoring mismatch with an N base
            bc = bc1
        else:
            mismatches += 1
            # Will resolve conflict if one of the barcodes is an *exact* match to a barcode on our list and the other is not.
            sample_1 = samples.get(bc1, 'unknown')
            sample_2 = samples.get(bc2, 'unknown')
            return  sample_1 if sample_1 == sample_2 else 'conflict' 
    if bc in samples:
        return samples[bc]
    else:
        matches = {sample for sample, matcher in sample_matchers.items() if matcher.find(bc) != -1}
        if len(matches) == 1:
            errors +=1 
            return matches.pop()
        else:
            return 'unknown'

extras = ['conflict', 'unknown']
filenames = samples.tolist()+extras

file_tallies = defaultdict(int)

fastq_ext = '.fastq{:}'.format('' if args.derep or args.compression == 'none' else '.'+args.compression)
out_files = {filename:[ smart_open(os.path.join(args.forward_read_dir, filename)+fastq_ext, 'wb', makedirs=True),
                        smart_open(os.path.join(args.reverse_read_dir, filename)+fastq_ext, 'wb', makedirs=True)] for filename in filenames}

with smart_open(args.forward_reads) as forward_file, smart_open(args.reverse_reads) as reverse_file:
    for i, (forward_line, reverse_line) in enumerate(zip(forward_file, reverse_file)):
        if i%4 == 0:
            filename = identify_barcode(forward_line.decode('ascii').split(':')[-1][:length], reverse_line.decode('ascii').split(':')[-1][:length])
            file_tallies[filename] += 1
            forward_out, reverse_out = out_files[filename]
        forward_out.write(forward_line)
        reverse_out.write(reverse_line)

for forward, reverse in out_files.values():
    forward.close(), reverse.close()

S = pd.Series(file_tallies)
total = S.sum()
successes = total - S[extras].sum()

Log("""Split {:,} reads ({:.2%} successfully).')
{:.2%} of forward-reverse barcodes matched perfectly.
{:.2%} of barcodes were saved by tolerating single-nucleotide errors in the barcode. 
Reads were split into the following files:
{:}""".format(
total, successes/total, 1 - mismatches/total, errors/total,
(S/total).to_string(float_format='{:.2%}'.format)
))

