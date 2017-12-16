#!/usr/bin/env python3

import argparse
import os
from subprocess import Popen, PIPE
import pandas as pd
from tuba_seq.fastq import fastqDF
from tuba_seq.shared import logPrint
from tuba_seq.pmap import CPUs

parser = argparse.ArgumentParser(description="Run PEAR (Illumina Paired-End reAd mergerR) on all forward and reverse fastq files.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("forward_read_dir", help="Directory containing forward read files.")
parser.add_argument("reverse_read_dir", help="Directory containing reverse read files (will be mated by sample name).")
parser.add_argument("-m", '--merge_dir', default='merged_reads', help="Directory to saved merged fastq files")
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multi-threaded operation')
parser.add_argument('-c', '--cmd', default='pear', help='Name of PEAR PATH/executable.')
parser.add_argument('-u', '--uncompressed', action='store_true', help='Avoid PHRED compression.')
parser.add_argument('-k', '--keep', action='store_true', help='Keep unassembled read files.')
parser.add_argument('-n', '--min_length', type=int, help="Minimum length for a merged sequence.", 
    default=(22*2)+29) #Default = minimum length necessary for successful preprocessing, under default conditions.  

# When merging paired-end reads, a combined PHRED score must be assigned to each
# nucleotide. From the trained error model, we can adjudicate whether the combined 
# PHRED score faithfully describes the true error rate. For example, a combined
# PHRED score of 40 *ought* to yield the expected base 99.99% of the time. PEAR 
# simply adds PHRED scores and appears to be too aggressive, i.e. a PHRED score 
# of 40 is correct <99.99% of the time. Thus, the PHRED scores are compressed:
# PHRED scores 1-60 are divided by three, while PHRED scores 
# above 60 are not compressed, but reduced by 40--to generate an uninterrupted 
# mapping. DADA2 is designed to automatically determine the true error rate of 
# each PHRED score, so none of this matters too much; however, DADA2 makes a 
# linear assumption during read-dereping that is inappropriate when the true 
# error rate and reported error rates disagree by several decades--as they will 
# if the PHRED scores aren't (crudely) compressed. 

ASCII_BASE = 33
attenuation_rate = 3
attenuation_cap = 60
max_attenuated = int(attenuation_cap/attenuation_rate)
max_PHRED = ASCII_BASE + 94
PHRED_compressor = {ASCII_BASE + i:(ASCII_BASE + int(i/attenuation_rate)) if i < attenuation_cap else (ASCII_BASE + i - attenuation_cap + max_attenuated) for i in range(max_PHRED)}

fastq_ext = '.fastq'

args = parser.parse_args()
Log = logPrint(args)

os.makedirs(args.merge_dir, exist_ok=True)

forward_files = {f for f in os.listdir(args.forward_read_dir) if fastq_ext in f}
reverse_files = {f for f in os.listdir(args.reverse_read_dir) if fastq_ext in f}

matches = forward_files & reverse_files
Log("Found {:} matching files.".format(len(matches)))

forward_only = forward_files - reverse_files
reverse_only = reverse_files - forward_files

if len(forward_only) + len(reverse_only) == 0:
    Log("Matched all fastq files")
else:
    if len(forward_only) > 0:
        Log("Could not find a reverse fastq file for:")
        for f in forward_only:
            Log(forward_only)
    if len(reverse_only) > 0:
        Log("Could not find a forward fastq file for:")
        for f in reverse_only:
            Log(forward_only)

stats = {'Assembled reads', 'Discarded reads', 'Not assembled reads'}
tallies = dict()

suffixes = {'assembled', 'discarded', 'unassembled.forward', 'unassembled.reverse'}
for File in matches:
    sample = File.split('.fastq')[0]
    output_file = os.path.join(args.merge_dir, sample)
    options = { '-f':os.path.join(args.forward_read_dir, File),
                '-r':os.path.join(args.reverse_read_dir, File),
                '-o':output_file,
                '-n':args.min_length,
                '-j':CPUs if args.parallel else 1,  # No. of threads to use
                '-c':max_PHRED}                             
    
    command = [args.cmd]+[str(s) for item in options.items() for s in item]
    Log('Analyzing {:} with command:\n{:}'.format(sample, ' '.join(command)))
    output = Popen(command, stdout=PIPE, stderr=PIPE).communicate()[0].decode('ascii')
    
    tallies[sample] = pd.Series({stat:int(line.partition(':')[2].partition('/')[0].replace(',', '')) 
                        for line in output.splitlines() for stat in stats if stat+' ...' in line}, name='Totals')
    
    if args.uncompressed:
        continue
    reads = fastqDF.from_file(output_file+".assembled.fastq", fake_header=False, use_Illumina_filter=False)
    reads['QC'] = reads['QC'].str.translate(PHRED_compressor)
    reads.write(output_file)
    if not args.keep:
        for suffix in suffixes:
            os.remove('{:}.{:}.fastq'.format(output_file, suffix))

tallies = pd.DataFrame(tallies).T
tallies.index.names = ['Sample']
Log(tallies.to_string())
Totals = tallies.sum()

Log('Summary of Read Merging', True, header=True)
df = pd.DataFrame(dict(Totals=Totals, Fraction=Totals/Totals.sum()))
Log(df.to_string(float_format='{:.3%}'.format), True)

