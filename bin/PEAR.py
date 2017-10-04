#!/usr/bin/env python3

#from params import barcode_length, alignment_flank
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

#PEAR_params =pd.Series({#'-v':['min_overlap', barcode_length, "Minimum overlap region length for a merged sequence."], 
                        #'-n':['min_length', barcode_length+2*alignment_flank, "Minimum length for a merged sequence."],
#                        '-c':['phred_cap', 0, 'Cap to PHRED score (0 = no cap).'],
#                        '-b':['base', 1, 'Base PHRED quality score.']}, index=['name', 'default', 'help'])

#PEAR_group = parser.add_argument_group('PEAR', 'Arguments passed directly to PEAR, as defined in its documentation.')
#for char, S in PEAR_params.items():
#    PEAR_group.add_argument(char, '--'+S['name'], **(S[['default', 'help']].to_dict()))

ASCII_BASE = 33
attenuation_rate = 3
attenuation_cap = 60
max_attenuated = int(attenuation_cap/attenuation_rate)
max_PHRED = 48

fastq_ext = '.fastq'

args = parser.parse_args()
Log = logPrint(args)

threads = CPUs if args.parallel else 1

os.makedirs(args.merged_dir, exist_ok=True)

forward_files = {f for f in os.listdir(parser.forward_read_dir) if fastq_ext in f}
reverse_files = {f for f in os.listdir(parser.reverse_read_dir) if fastq_ext in f}

matches = forward_files & reverse_files
Log("Found "+len(matches)+" matching files.")

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

PHRED_compressor = {ASCII_BASE + i:(ASCII_BASE + int(i/attenuation_rate)) if i < attenuation_cap else (ASCII_BASE + i - attenuation_cap + max_attenuated) for i in range(max_PHRED)}
Totals = pd.Series({'Assembled reads':0, 'Discarded reads':0, 'Not assembled reads':0})

for File in matches:
    output_file = os.path.join(args.merge_dir, File)
    options = { '-f':os.path.join(args.forward_read_dir, File),
                '-r':os.path.join(args.reverse_read_dir, File),
                '-o':output_file,
                '-j':threads,
                '-c':0,
                '-b':1}
    
    #options.update(PEAR_params['name'].apply(args.__dict__.get).to_dict())

    command = [args.PEAR_CMD]+[str(s) for item in options.items() for s in item]
    Log('Analyzing '+File+' with command:')
    Log(' '.join(command))
    output = Popen(command, stdout=PIPE, stderr=PIPE).communicate()[0].decode('ascii')
    
    for line in output.splitlines():
        for stat in Totals.keys():
            if stat in line:
                Log(line)
                Totals[stat] += int(line.partition(':')[2].partition('/')[0].replace(',', '')) 
    
    if args.uncompressed:
        continue
    reads = fastqDF.from_file(output_file, fake_header=False, use_Illumina_filter=False)
    reads['QC'] = reads['QC'].str.translate(PHRED_compressor)
    reads.write(output_file)

df = pd.DataFrame(dict(Totals=Totals, Fraction=Totals/Totals.sum()))
Log(df.to_string(float_format='{.3%}'.format), True)

