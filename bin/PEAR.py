#!/usr/bin/env python3

import argparse
import os
from subprocess import Popen, PIPE
import pandas as pd
from tuba_seq.shared import logPrint
from tuba_seq.pmap import CPUs

parser = argparse.ArgumentParser(description="Run PEAR (Illumina Paired-End reAd mergerR) on all forward and reverse fastq files.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("forward_reads", help="Directory containing forward read files.")
parser.add_argument("reverse_reads", help="Directory containing reverse read files (will be mated by sample name).")
parser.add_argument("-m", '--merge_dir', default='merged_reads', help="Directory to saved merged fastq files")
parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-p', '--parallel', action='store_true', help='Multi-threaded operation')
parser.add_argument('-c', '--cmd', default='pear', help='Name of PEAR PATH/executable.')
parser.add_argument('-k', '--keep', action='store_true', help='Keep unassembled read files.')
parser.add_argument('-n', '--min_length', type=int, help="Minimum length for a merged sequence.", 
    default=(22*2)+29) #Default = minimum length necessary for successful preprocessing, under default conditions.  
parser.add_argument('--memory', type=str, default='16G', help='Memory to use. K,M,G are possible suffixes.')
parser.add_argument('--max_PHRED', type=int, default=72, help='Maximum PHRED score allowed.')

fastq_ext = '.fastq'

args = parser.parse_args()
Log = logPrint(args)

os.makedirs(args.merge_dir, exist_ok=True)

single_file = '.fastq' in args.forward_reads
if single_file:
    assert '.fastq' in args.reverse_reads, """First argument is a single FASTQ file; however, the second argument isn't a FASTQ file. 
You must either give single FASTQ files to this program or entire directories."""

if not single_file:
    forward_files = {f for f in os.listdir(args.forward_reads) if fastq_ext in f}
    reverse_files = {f for f in os.listdir(args.reverse_reads) if fastq_ext in f}

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
else:
    args.forward_reads, File = os.path.split(args.forward_reads)
    args.reverse_reads, reverse_file = os.path.split(args.reverse_reads)
    assert reverse_file == File, 'Forward & Reverse files must have the same basename.'
    matches = [File]

stats = {'Assembled reads', 'Discarded reads', 'Not assembled reads'}
tallies = dict()

suffixes = {'discarded', 'unassembled.forward', 'unassembled.reverse'} 

for File in matches:
    sample = File.split(fastq_ext)[0]
    output_file = os.path.join(args.merge_dir, sample)
    options = { '-f':os.path.join(args.forward_reads, File),
                '-r':os.path.join(args.reverse_reads, File),
                '-o':output_file,
                '-n':args.min_length,
                '-j':CPUs if args.parallel else 1,  # No. of threads to use
                '-c':args.max_PHRED,
                '-y':args.memory}                             
    
    command = [args.cmd]+[str(s) for item in options.items() for s in item]
    Log('Analyzing {:} with command:\n{:}'.format(sample, ' '.join(command)))
    output = Popen(command, stdout=PIPE, stderr=PIPE).communicate()[0].decode('ascii')
    tallies[sample] = pd.Series({stat:int(line.partition(':')[2].partition('/')[0].replace(',', '')) 
                        for line in output.splitlines() for stat in stats if stat+' ...' in line}, name='Totals')
    
    os.rename(output_file+'.assembled.fastq', output_file+'.fastq')
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

