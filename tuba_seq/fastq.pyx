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

NW_kwargs = dict(match=2, mismatch=1, gap_open=3, gap_extend=1)
from striped_smith_waterman import SW
sw = SW(**NW_kwargs)

def cprint(s): print(s.decode('ascii'))

from shared import smart_open
class IterFASTQ(object): 
    def __iter__(self): return self
    
    def __init__(self, filename): 
        self.filename = filename
        self.f = smart_open(self.filename)

    def __next__(self):
        header = self.f.readline()
        if not header:
            self.f.close()
            raise StopIteration
        dna = self.f.readline()
        self.f.readline()
        QC = self.f.readline()
        if not QC:
            self.f.close()
            raise RuntimeError("Input FASTQ file was not 4x lines long")
        return header, dna, QC

cdef:
    bytes LINE_3 = b"\n+\n"
    bytes END = b'\n'
    bytes ILLUMINA_FAILED_FILTER = b':Y:'
    bytes c_gap = b'-'
    bytes c_N = b'N'
    int o_gap = ord(c_gap)
    int o_N = ord(c_N)
    char char_gap = c_gap[0]
    char char_N = c_N[0]
    int DNA_to_int[256]

for i, nuc in enumerate(b'ACGTN-'):
    DNA_to_int[nuc] = i


class MasterRead(object):
    possible_outcomes = ['Filtered', 'Unaligned', 'Wrong Barcode Length', 'Residual N', 'Insufficient Flank', 'Clustered']
    MAX_READ_LENGTH = 300

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
        full_length = len(self.full)

        if args.trim == 'symmetric':
            self.pre_slice = slice(start - aft_length, full_length) if start > aft_length else slice(0, full_length - (aft_length - start))
        elif ',' in args.trim:
            pos = args.trim.split(',')
            if len(pos) != 2:
                raise RuntimeError("--trim must assign a start and stop position")
            self.pre_slice(int(pos[0]), int(pos[1]))
        else:
            self.pre_slice = slice(0, full_length)
       
        self.barcode_length = stop - start
        
        self.c_ref = self.ref.encode('ascii')
        self.c_train = self.c_ref[:self.training_flank] + self.c_ref[-self.training_flank:]
        self.max_score = sw.char_score(self.c_ref, self.c_ref)
        self.min_align_score = args.min_align_score
        self.min_int_score = int(np.ceil(args.min_align_score*self.max_score))
        
    def iter_fastq(self, input_fastq_iter, filenames):
        scores = pd.Series(np.zeros(self.max_score+1, dtype=int), index=pd.Index(np.linspace(0,1,num=self.max_score+1), name='Score'), name='Occurrences')
        bad_barcode_lengths = pd.Series(np.zeros(self.MAX_READ_LENGTH, dtype=int), index=pd.Index(np.arange(self.MAX_READ_LENGTH), name='Length'), name='Occurrences')
        cdef:
            int BL = self.barcode_length
            int TF = self.training_flank
            int CF = self.cluster_flank
            int start, stop
            int score
        
            int Filtered = 0
            int Residual_N = 0
            int Insufficient_Flank = 0
            int Clustered = 0
            long [:] score_view = scores.values
            long [:] bc_length_view = bad_barcode_lengths.values
            int qc_i
            int unaligned_counter = 0

        with smart_open(filenames[0], 'wb', makedirs=True) as training_file, smart_open(filenames[1], 'wb', makedirs=True) as cluster_file, smart_open(filenames[2], 'wb', makedirs=True) as unaligned_file:
            for header, DNA, QC in input_fastq_iter: 
                DNA = DNA[self.pre_slice]
                QC = QC[self.pre_slice]
                if not QC:
                    raise RuntimeError("Input FASTQ file was not 4x lines long")
                if ILLUMINA_FAILED_FILTER in header:
                    Filtered += 1
                    continue
                start, stop = sw.find_N_start_stop(DNA, self.c_ref)
                DNA_scoring = DNA[start - self.alignment_flank:stop + self.alignment_flank]
                score = sw.char_score(DNA_scoring, self.c_ref)
                score_view[score] += 1
                if score < self.min_int_score:
                    unaligned_counter += 1
                    unaligned_file.write(DNA+END)
                    continue
                if abs((stop - start) - BL) > self.allowable_deviation:
                    bc_length_view[stop - start] += 1
                    continue
                T_s1 = slice(start - TF, start)
                T_s2 = slice(stop, stop+TF)
                training_DNA = DNA[T_s1]+DNA[T_s2]
                if len(training_DNA) == 2*TF and c_N not in training_DNA:
                    tQC = QC[T_s1]+QC[T_s2]
                    assert len(tQC) == len(training_DNA), '{:} {:} {:} {:}'.format(len(tQC), len(training_DNA), tQC, training_DNA)
                    training_file.write(header+training_DNA+LINE_3+tQC+END)
        
                if c_N in DNA:
                    DNA = sw.fill_Ns(DNA, self.c_ref)
                cluster_DNA = DNA[start - CF:start+BL+CF]
                if c_N in cluster_DNA:
                    Residual_N += 1
                    continue
                if len(cluster_DNA) != (BL+ 2*CF):
                    Insufficient_Flank += 1
                else:
                    Clustered += 1
                    cQC = QC[start-CF:start+BL+CF]
                    cluster_file.write(header+cluster_DNA+LINE_3+cQC+END)
        statistics = pd.Series([Filtered,   scores.iloc[:self.min_int_score].sum(),   bad_barcode_lengths.sum(),   Residual_N,   Insufficient_Flank,   Clustered], 
                            index=pd.Index(self.possible_outcomes))
        return statistics, scores, bad_barcode_lengths

import regex as re
class Mismatcher(object):
    def __init__(self, substring, n_substitutions=1):
        self.pattern_obj = re.compile("("+substring+'){s<='+str(n_substitutions)+'}')

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

def get_QC_map():
    return np.array(QC_map)[:128]

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

