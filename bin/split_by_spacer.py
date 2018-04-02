#!/usr/bin/env python3

import argparse
from tuba_seq.fastq import Mismatcher, IterFASTQ
from tuba_seq.shared import smart_open, logPrint
import pandas as pd
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser(description="Split paired-end read files by Illumina indecies.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('spacer_file', type=str, 
    help='''CSV file with the following columns (other columns OK).: 
`input_file`, `forward_spacer_length`, `reverse_spacer_length`, `sample`
''')
parser.add_argument("--input", type=Path, default=Path('merged_reads'),
    help='FASTQ file or directory of FASTQ files')
parser.add_argument("--output_dir", type=Path, default=Path('split_reads'),
    help='Directory for output FASTQ files.')
parser.add_argument("--failed_dir", type=Path, default=Path("unsplit_reads"),
    help='Directory for reads that could not be split')
parser.add_argument('--forward', type=str, default='ATGCGCACGT', 
    help='Forward string')
parser.add_argument('--reverse', type=str, default='ATGTCCAAGA'[::-1],
    help='Reverse string (read in the forward direction)')
parser.add_argument('-p', '--parallel', action='store_true', 
    help='Multithreaded operation')
parser.add_argument('--compression', default='gz', choices=['bz2', 'gz', 'lzma', 'none'], 
    help='Compression algorithm for output.')
parser.add_argument('--indel', default=1, type=int, 
    help='Tolerable deviation from expected spacer.')
parser.add_argument('--substitutions', default=3, type=int, 
    help='Tolerable substitutions from expected string')
parser.add_argument('--use_reverse', action='store_true', 
    help='Use the reverse N-length to more-conservatively assign reads.')
parser.add_argument('--reverse_inset', type=int, default=4, help='Inset of reverse spacer from end of the read (which usually includes the 5nt T7 index)')
###############################################################################


args = parser.parse_args()
Log = logPrint(args)

samples = pd.read_csv(args.spacer_file).set_index(['input_file', 'forward_spacer_length', 'reverse_spacer_length'])['sample']

fastq_ext = '.fastq'
if args.compression != 'none':
    fastq_ext += '.'+args.compression
if args.parallel:
    from tuba_seq.pmap import pmap as map

from collections import defaultdict
class StringFinder(object):
    def __init__(self, string, acceptable_spacers): 
        acceptable_spacers = np.array(list(acceptable_spacers))
        acceptable_spacers.sort()
        if len(acceptable_spacers) > 1:
            min_diff = abs(np.diff(acceptable_spacers)).min() 
            assert min_diff > 2*args.indel, "The tolerated indel is {:} nts, yet there are spacers that differ by only {:} nts".format(args.indel, min_diff)
        
        self.spacer_lengths_dict = {acceptable:spacer for spacer in acceptable_spacers for acceptable in range(spacer-args.indel, spacer+args.indel+1) }
        self.bad_lengths = defaultdict(int)
        self.successes = 0
        self.finder = Mismatcher(string, n_substitutions=args.substitutions)
        self.acceptable_spacers = acceptable_spacers

    def find(self, DNA_seq):
        start = self.finder.find(DNA_seq)
        if start in self.spacer_lengths_dict:
            self.successes += 1
            return self.spacer_lengths_dict[start]
        self.bad_lengths[start] += 1
        return 'Failed'

    def summarize(self):
        bad_lengths = self.bad_lengths.copy()
        no_match = bad_lengths.pop(-1)
        worst_length = max(bad_lengths, key=bad_lengths.get) 
        return pd.Series({  'No Match':no_match, 
                            'Bad Length':sum(bad_lengths.values()), 
                            'Success':self.successes,
                            'Most Common Bad Length ({:})'.format(worst_length):self.bad_lengths[worst_length]})

N_bad_combos = 3
Files = args.input.glob('*.fastq*') if args.input.is_dir() else [args.input]
def split_fastq(filename):
    bad_combos = defaultdict(int)
    reads = 0
    input_name = filename.stem.split('.')[0]
    spacers = samples.loc[filename.name].reset_index()
    
    mapper = spacers.set_index(['forward_spacer_length', 'reverse_spacer_length'])['sample']
    forward_spacers = frozenset(spacers['forward_spacer_length'].values)
    reverse_spacers = frozenset(spacers['reverse_spacer_length'].values)
    forward_finder = StringFinder(args.forward, forward_spacers)
    reverse_finder = StringFinder(args.reverse[::-1], reverse_spacers)
    output_files = {sample:smart_open(args.output_dir / (sample + fastq_ext), 'wb', makedirs=True) for sample in spacers['sample']}
    output_files['Failed'] = smart_open(args.failed_dir / (input_name + fastq_ext), 'wb', makedirs=True)
    destinations = dict(zip(output_files.keys(), len(output_files)*[0]))
    for bheader, bDNA, bQC in IterFASTQ(filename):
        DNA = bDNA.decode('ascii')
        spacers = forward_finder.find(DNA), reverse_finder.find(DNA[::-1]) 
        if spacers in mapper:
            fastq = mapper[spacers]
        elif spacers[0] in forward_spacers:
            if spacers[1] in reverse_spacers:
                bad_combos[spacers] += 1
            fastq = 'Failed' if args.use_reverse else mapper[spacers[0]].iloc[0] 
        else: 
            fastq = 'Failed' 
        destinations[fastq] += 1 
        output_files[fastq].write(bheader+bDNA+b'+\n'+bQC)
        reads += 1
        if reads >= 100000:
            break

    for output_file in output_files.values(): 
        output_file.close()
    
    failures = destinations.pop("Failed")
    successes = reads - failures
    bad_combos = pd.Series(bad_combos)
    outcomes = pd.concat({  'Forward' : forward_finder.summarize(),
                            'Reverse' : reverse_finder.summarize(),
                            'Overall' : pd.Series({
                                'Bad Combo'    : int(bad_combos.sum()), 
                                'Success'      : successes,
                                'reads'        : reads})
                        })
    percents = outcomes/reads    
    verbose = percents['Overall', 'Success'] < 0.9
    Log("Successfully split {:.2%} of {:,} reads in {:} into:".format(percents['Overall', 'Success'], reads, input_name), True)
    reverse_mapper = mapper.reset_index().set_index("sample")
    for name, destination in destinations.items():
        Log('{:>20}: {:.1%} of successful reads, {:}N/{:}N spacers.'.format(name, destination/successes, *reverse_mapper.loc[name]), True)
    Log(percents.to_string(float_format='{:.2%}'.format), verbose)
    #TopN_bad_combos = (bad_combos.sort_values(ascending=False)/bad_combos.sum()).iloc[:N_bad_combos]
    #Log('Top {:} bad combos (% of total bad combos):\n'.format(N_bad_combos)+TopN_bad_combos.to_string(float_format='{:.2%}'.format), verbose)
    for condition, stat in outcomes.keys():
        if 'Most Common Bad Length' in stat:
            outcomes.pop((condition, stat))
    return input_name, outcomes.sort_index()

outcomes = pd.DataFrame(dict(map(split_fastq, Files))).T
reads = outcomes.pop(('Overall', 'reads')).sum()
totals = outcomes.sum()
T_percents = totals/reads

Log("Successfully split {:.2%} of {:,} reads in {:} files.".format(T_percents['Overall', 'Success'], reads, len(outcomes)), True)
Log(T_percents.to_string(float_format='{:.3%}'.format), True)

