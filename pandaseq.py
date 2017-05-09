#!/bin/bash
from params import panda_seq
import argparse
import os

parser = argparse.ArgumentParser(description="Run pandaseq on all forward and reverse fastq files.")
forward_file_glob = "raw/stacked/forward/{sample}.fastq*"
reverse_file_glob = "raw/stacked/reverse/{sample}.fastq*"

parser.add_argument("-v", "--verbose", help='Output more Info', action="store_true")
parser.add_argument('-L', '--log', action='store_true', help='Log output')
args = parser.parse_args()

#algorithm = 'stitch'        
algorithm = 'simple_bayesian'
os.makedirs(panda_seq['merge_dir'], exist_ok=True)

def collect_sample_dict(file_glob, glob_chars='?*\[\]'):
    from glob import glob
    import re
    files = glob(file_glob.format(sample='*'))
    extract_str = [s for s in re.split('['+glob_chars+']', file_glob) if '{sample}' in s][0]
    out_dict = {}
    for f in files:
        finds = re.findall(extract_str.format(sample='(.+)'), f)
        if len(finds) != 1:
            print("Could not extract sample name from file:", f)
            raise ValueError
        sample = finds[0]
        if sample in out_dict:
            print("Found multiple files with same sample name:", f, 'and', out_dict[sample])
            raise RuntimeError
        out_dict[sample] = f
    return out_dict

forward_files = collect_sample_dict(forward_file_glob)
reverse_files = collect_sample_dict(reverse_file_glob)

if len(forward_files) == len(reverse_files):
    print("Found", len(forward_files), "samples.")
else:
    print("Found", len(forward_files), "forward files and", len(reverse_files), "reverse files.")

map_params_to_pandaseq_option = dict(
                    l='min_length',       # Minimum length for a sequence.
                    t='threshold',        # The minimum probability that a sequence must have to assemble and, if used, match a primer.
                    O='min_overlap')      # Minimum overlap region length for a sequence. (0 to use read length.)

extra_params = {s:panda_seq[param_name] for s, param_name in map_params_to_pandaseq_option.items()}

def merge_files(sample, forward_file, reverse_file, cmd='pandaseq'):
    from tuba_seq.pmap import CPUs
    from subprocess import Popen, PIPE
    options = dict( f=forward_file,
                    r=reverse_file,
                    W=panda_seq['merge_dir']+sample+'.fastq.bz2',
                    T=CPUs,
                    A=algorithm)# algorithm:parameters	Select the algorithm to use for overlap selection and scoring.
    options.update(extra_params)
    if args.log:
        options['G'] = sample+'.LOG.bz2'

    str_options = {'-'+k:str(v) for k, v in options.items()}
    return Popen([cmd, '-a', '-F']+[s for item in str_options.items() for s in item], stdout=PIPE, stderr=PIPE).communicate()

for sample, forward_file in forward_files.items():
    if sample in reverse_files:
        output = merge_files(sample, forward_file, reverse_files[sample]) 
        #if not args.log:
            #statlines = output[1].splitlines()[-8:-1]
            #stats = {name:int(count) for __, __, name, count in map(lambda b: b.decode('ascii').split('\t'), statlines)} 
            #print("Merged {:.1%} of {:.1}M reads in sample {:}.".format(stats['OK']/stats['READS'], stats['READS']*1e-6, sample))
    else:
        print("Could not find a reverse FASTQ for", forward_file)
