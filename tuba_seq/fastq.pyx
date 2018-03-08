"""Low-level, efficient functionality to handle FASTQ files.

There are three classes within this module:

1) fastqDF (subclass of pandas.DataFrame)
    Reads, processes (slices, queries, etc), and writes FASTQ files.

2) MasterRead
    Can identify a `MasterRead` from amplicon pileups and then align/score reads
    against this `MasterRead`. 

3) singleMismatcher
    Identifies single-mismatch-tolerant substrings within a sequence.

"""

import pandas as pd
import os
from collections import OrderedDict
import numpy as np
cimport numpy as np

class fastqDF(pd.DataFrame):
    essential_columns = ['header', 'DNA', 'QC']
    def __init__(self, *args, **kargs):
        super(fastqDF, self).__init__(*args, **kargs)

    def select_reads(self, ix):
            reduced = fastqDF(self.loc[ix, :])
            reduced.info = self.info
            return reduced

    @classmethod
    def from_file(cls, filename, fake_header=True, use_Illumina_filter=True, check_lengths=True):
        self = cls(pd.DataFrame(
            data=pd.read_table(filename, names=['__'], encoding='ascii', quoting=3, dtype=str)['__'].values.reshape((-1, 4))[:, [True, True, False, True]], 
            columns=cls.essential_columns))
        if use_Illumina_filter:
            self.select_reads(self['header'].str.contains(":N:"))       # Y = failed filter, N = passed filter
        sample_header = self['header'].iloc[0]
        split = sample_header.split()
        self.info = pd.Series(dict(zip(['Instrument', 'Run', 'Flowcell', 'Lane'], split[0].split(':'))))
        if len(split) == 2:
            h2 = split[1]
            self.info['Paired-end'] = h2.split(':')[0] == '1'
        self.info['Index'] = sample_header.split(':')[-1]
        self.info['Initial Reads'] = len(self)
        self.info['Filtered'] = use_Illumina_filter
        if fake_header:
            self.info['Fake Header'] = sample_header
            self.pop('header') 
        self.info['Sample'] = os.path.basename(filename.partition('.fastq')[0])
        if check_lengths:
            DNA_L = self['DNA'].str.len()
            QC_L = self['QC'].str.len()
            matched = DNA_L == QC_L
            if not matched.all():
                unmatched = self.loc[~matched]
                print('In FASTQ file', filename, 'the following reads exhibit DNA sequences and Phred Scores with inconsistent lengths:') 
                print(unmatched)
                raise RuntimeError("Inconsistent lengths in FASTQ File")
        return self

    def query(self, s, **kargs):
        return super(fastqDF, self).query(s, **kargs)

    def isDegenerate(self):
        return self['DNA'].str.contains('N') 

    def drop_abnormal_lengths(self):
        lengths = self['DNA'].str.len()
        keep = lengths == int(lengths[:999].median())
        return self.select_reads(keep)
    
    def co_slice(self, *args, **kargs):
        if args:
            kargs['start'] = args[0].start
            kargs['stop'] = args[0].stop
        self['DNA'] = self['DNA'].str.slice(**kargs)
        self['QC'] = self['QC'].str.slice(**kargs)
        return self

    def vector_slice(self, slices, second_slice=None):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            initial_seqs = self[['DNA', 'QC']].values

            self.loc[:, ['DNA', 'QC']] = np.array([[dna[Slice], qc[Slice]] for (dna, qc), Slice in zip(initial_seqs, slices)])

            if second_slice is not None:
                second_slice = np.array([[dna[Slice], qc[Slice]] for (dna, qc), Slice in zip(initial_seqs, slices)])
                self['DNA'] = self['DNA'].str.cat(second_slice[:, 0]) 
                self['QC' ] = self['QC' ].str.cat(second_slice[:, 1])
        return self

    def expected_errors(self):
        return self['QC'].apply(_expected_errors)

    def write(self, filename):
        if 'Fake Header' in self.info:
            self.insert(1, 'header', self.info['Fake Header'])
        if '__plus_sign__' not in self.columns:
            self.insert(2, '__plus_sign__', len(self)*['+'])
        filename = filename.split('.fastq')[0] + '.fastq'
        self.to_csv(filename + '.gz', 
                    columns=self.essential_columns[0:2] + ['__plus_sign__'] + self.essential_columns[2:3], 
                    compression='gzip',
                    header=False,
                    index=False,
                    sep='\n')

    def construct_read_set(self, dna, QCs):
        s = self.info['Fake Header']+'\n'+dna+'\n+\n'
        return (s+('\n'+s).join(QCs)+'\n').encode('ascii')

NW_kwargs = dict(match=6, mismatch=-3, gap_open=-12, gap_extend=-1)

from seq_align import NW
nw = NW(**NW_kwargs)
nw.add_neutral_N()

def cprint(s): print(s.decode('ascii'))

cdef bytes c_gap = '-'.encode('ascii') 
cdef bytes c_N = 'N'.encode('ascii') 
cdef bytes c_2 = '2'.encode("ascii")

from shared import smart_open
def iter_fastq(filename):
    with smart_open(filename) as f:
        while True:
            header = f.readline()
            if not header:
                break
            dna = f.readline()
            f.readline()
            QC = f.readline()
            if not QC:
                raise RuntimeError("Input FASTQ file was not 4x lines long")
            yield header, dna, QC

class MasterRead(object):
    def __init__(self, master_read, args):
        self.alignment_flank = args.alignment_flank
        self.training_flank = args.training_flank
        self.cluster_flank = args.cluster_flank if hasattr(args, 'cluster_flank') else args.allowable_deviation
        self.allowable_deviation = args.allowable_deviation

        self.full = master_read.replace('.', 'N')
        
        start = self.full.index("N")
        stop = self.full.rindex('N') + 1
        self.ref = self.full[start - self.alignment_flank:stop+self.alignment_flank]
        
        aft_length = len(self.full) - stop
        
        self.trim = min(args.trim, start - self.alignment_flank, aft_length - self.alignment_flank)
        if self.trim < 0:
            raise RuntimeError("The master read's flanks are shorter than `alignment_flank`")
            
        if args.symmetric_flanks:
            self.pre_slice = slice(self.trim + start - aft_length, len(self.full) - self.trim) if start > aft_length else slice(self.trim, len(self.full) - self.trim - aft_length + start)
        else:
            self.pre_slice = slice(self.trim, len(self.full) - self.trim) 
       
        self.barcode_length = stop - start
        
        if start - self.alignment_flank < self.pre_slice.start:
            raise ValueError("Insufficient forward flank ({:} nt) on Master Read to accommodate `alignment_flank` ({:} nt) & `trim` ({:} nt).".format(start, self.alignment_flank, self.trim))
        if stop + self.alignment_flank > self.pre_slice.stop:
            raise ValueError("Insufficient aft flank ({:} nt) on Master Read to accommodate `alignment_flank` ({:} nt) & `trim` ({:} nt).".format(aft_length, self.alignment_flank, self.trim))
        
        # MAYBE DELETE SOMETIME!
        assert len(self.full[self.pre_slice]) >= 2*self.alignment_flank + self.barcode_length 
        if args.symmetric_flanks:
            assert start - self.pre_slice.start == self.pre_slice.stop - stop, 'symmetric_flanks failed.' 
        ######

        self.c_ref = self.ref.encode('ascii')
        self.c_train = self.c_ref[:self.training_flank] + self.c_ref[-self.training_flank:]
        self.max_score = nw.char_score(self.c_ref, self.c_ref)
        self.min_align_score = args.min_align_score
        
    def find_start_stop(self, seq): 
        """score=score of [match, mismatch, gap_start, gap_extend]"""
        seq_align, ref_align, _score = nw.char_align(seq, self.c_ref) 
        cdef:
            int ref_start = ref_align.index(c_N)
            int ref_stop = ref_align.rindex(c_N) + 1
            int head_gaps = seq_align[0:ref_start].count(c_gap)
            int barcode_gaps = seq_align[ref_start:ref_stop].count(c_gap)
        return ref_start - head_gaps, ref_stop - head_gaps - barcode_gaps

    def repair_N(self, c_seq):
        seq_align, ref_align, _score = nw.char_align(c_seq, self.c_ref)
        return ''.encode('ascii').join([s if (s != c_N or r == c_gap) else r for s, r in zip(seq_align, ref_align) if s != c_gap])
    
    def tally_mutations(self, seq, QC):
        """Keep track of all the deviations from the reference and their QC score."""
        #assert len(seq) == len(QC)
        seq_align, ref_align, _score = nw.char_align(seq, self.c_train)
        qc_i = 0
        for s, r in zip(seq_align.decode('ascii'), ref_align.decode('ascii')):
            if s != '-':
                if s != 'N' and r != 'N' and r != '-':
                    ix = r+'2'+s, QC[qc_i] - 32
                    self.PHRED_tally[ix] = self.PHRED_tally.get(ix, 0) + 1
                qc_i += 1
    
    def iter_fastq(self, sample, filenames, input_fastq):
        from collections import defaultdict
        scores = defaultdict(int)
        cdef:
            int BL = self.barcode_length
            int TF = self.training_flank
            int CF = self.cluster_flank
            int start, stop
            double score
        
            int Filtered = 0
            int Unaligned = 0
            int Wrong_Barcode_Length = 0
            int Residual_N = 0
            int Insufficient_Flank = 0
            int Clustered = 0
            
            bytes LINE_3 = '\n+\n'.encode('ascii')
            bytes END = '\n'.encode('ascii')
            bytes ILLUMINA_FAILED_FILTER = ':Y:'.encode('ascii') 
              
        with smart_open(filenames[0], 'wb', makedirs=True) as training_file, smart_open(filenames[1], 'wb', makedirs=True) as cluster_file, smart_open(input_fastq) as input_file:
            while True:
                header = input_file.readline()
                if not header:
                    break
                DNA = input_file.readline()[self.pre_slice]
                input_file.readline()
                QC = input_file.readline()[self.pre_slice]
                if not QC:
                    raise RuntimeError("Input FASTQ file was not 4x lines long")
                if ILLUMINA_FAILED_FILTER in header:
                    Filtered += 1
                    continue
                if DNA in self.alignments:
                    start, stop, score = self.alignments[DNA]
                else:
                    start, stop = self.find_start_stop(DNA)
                    DNA_scoring = DNA[start - self.alignment_flank:stop + self.alignment_flank]
                    score = nw.char_score(DNA_scoring, self.c_ref)/self.max_score
                    if len(self.alignments) < self.max_alignments:
                        self.alignments[DNA] = start, stop, score
                scores[score] += 1
                if score < self.min_align_score:
                    Unaligned += 1
                    if hasattr(self, 'unaligned'):
                        self.unaligned[DNA] = 1 + self.unaligned.get(DNA, 0)
                    continue
                if abs((stop - start) - BL) > self.allowable_deviation:
                    Wrong_Barcode_Length += 1
                    continue
                T_s1 = slice(start - TF, start)
                T_s2 = slice(stop, stop+TF)
                training_DNA = DNA[T_s1]+DNA[T_s2]
                if c_N not in training_DNA and len(training_DNA) == 2*TF:
                    tQC = QC[T_s1]+QC[T_s2]
                    #assert len(tQC) == len(training_DNA), '{:}\n{:}'.format(training_DNA, tQC)
                    self.tally_mutations(training_DNA, tQC)                                     # To get my own tally
                    training_file.write(header+training_DNA+LINE_3+tQC+END)
                #if len(training_DNA) != 2*TF:
                #    print(2*TF - len(training_DNA))
                #if "N" in training_DNA:
                #    print('N')
                if c_N in DNA:
                    DNA = self.repair_N(DNA)
                cluster_DNA = DNA[start - CF:start+BL+CF]
                if c_N in cluster_DNA:
                    Residual_N += 1
                    continue
                if len(cluster_DNA) != (BL+ 2*CF):
                    Insufficient_Flank += 1
                else:
                    Clustered += 1
                    cQC = QC[start-CF:start+BL+CF]
                    assert len(cQC) == len(cluster_DNA), '{:}\n{:}'.format(cluster_DNA, cQC)
                    cluster_file.write(header+cluster_DNA+LINE_3+cQC+END)

        self.instruments[sample] = header.decode('ascii').split(':')[0]
        self.scores[sample] = pd.Series(scores)
        return pd.Series([Filtered,   Unaligned,   Wrong_Barcode_Length,   Residual_N,   Insufficient_Flank,   Clustered], 
         index=pd.Index(['Filtered', 'Unaligned', 'Wrong Barcode Length', 'Residual N', 'Insufficient Flank', 'Clustered']))

import regex as re
class singleMismatcher(object):
    def __init__(self, substring):
        self.pattern_obj = re.compile(substring+'{s<=1}')

    def find(self, searchstring):
        search_obj = self.pattern_obj.search(searchstring)
        return search_obj.start() if search_obj is not None else -1

def hamming_distance(a, b):
    from scipy.spatial.distance import hamming
    N = len(a)
    return round(N*hamming(np.fromiter(a, 'S1', N), np.fromiter(b, 'S1', N)))

class robustMatcher(object):
    def __init__(self, match_dict, max_errors=1):
        self.exact = 0
        self.recovered = 0
        self.conflict = 0
        self.unknown = 0
        self.match_dict = match_dict
        self._recover_str = ' '.joinmatch_dict.values()
        self._ext_str = '{{s<={:}}}'.format(max_errors)

    def match(self, s):
        if s in self.match_dict:
            self.exact += 1
            return self.match_dict[s]
        recoveries = re.findall(s+self._ext_str, self._recover_str)
        L = len(recoveries)

        if L == 1:
            self.recovered += 1
            return self.match_dict[recoveries[0]]
        if L == 0:
            self.unknown += 1
            return "unknown"
        self.conflict += 1
        return "conflict"

    def summary(self):
        return pd.Series({outcome:getattr(self, outcome) for outcome in ['exact', 'recovered', 'unknown', 'conflict']})

cdef double QC_map_C[256] 
cdef double [:] QC_map = QC_map_C
cdef double dp = 0.7943282347
QC_map[33] = 1.0
for i in range(34, 128):
    QC_map[i] = dp*QC_map[i - 1]

def _expected_errors(QC):
    Q_score_s = QC.encode('ascii')
    cdef char *Q_scores = Q_score_s
    cdef int L = len(QC)
    cdef int q_score, i
    cdef double errors = 0
    for i in range(L):
        q_score = <int>Q_scores[i]
        if q_score < 33 or q_score > 75:
            raise ValueError("Q score must be between 33 and 75. It was {:} (i.e. {:})".format(q_score, QC.encode('ascii')[i]))
            
        errors += QC_map[q_score]
    return errors

