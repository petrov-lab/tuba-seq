import pandas as pd
from collections import OrderedDict
import numpy as np
cimport numpy as np

class fastqDF(pd.DataFrame):
    essential_columns = ['header', 'DNA', 'QC']
    def __init__(self, filename, immediate_trunc):
        linear = pd.read_table(filename, names=['__'])['__'] 
        header = linear.loc[::4].values
        DNA = linear.loc[1::4].str.slice(immediate_trunc.start, immediate_trunc.stop).values
        QC  = linear.loc[3::4].str.slice(immediate_trunc.start, immediate_trunc.stop).values
        pd.DataFrame.__init__(self, OrderedDict(zip(self.essential_columns, [header, DNA, QC])))
        init_reads = len(self)
        keep = self['header'].str.contains(":N:")   # Y = filtered, N = passed filter
        self = self[keep]
        self.short_filename = filename.partition('.fastq')[0].split('/')[-1]

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def query(self, s, **kargs):
        self = super(fastqDF, self).query(s, **kargs)
        return self

    def isDegenerate(self):
        return self['DNA'].str.contains('N') 

    def dropDegen(self):
        droped = self.loc[-self.isDegenerate(), :]
        droped.__class__ = fastqDF
        return droped

    def drop_abnormal_lengths(self):
        lengths = self['DNA'].str.len()
        keep = lengths == int(lengths.median())
        droped = self.loc[keep, :]
        droped.__class__ = fastqDF
        return droped
    
    def vslice(self, slices, second_slice=None, enforce=True):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            initial_seqs = self[['DNA', 'QC']].values

            self.loc[:, ['DNA', 'QC']] = np.array([[dna[Slice], qc[Slice]] for (dna, qc), Slice in zip(initial_seqs, slices)])

            if second_slice is not None:
                second_slice = np.array([[dna[Slice], qc[Slice]] for (dna, qc), Slice in zip(initial_seqs, slices)])
                self['DNA'] = self['DNA'].str.cat(second_slice[:, 0]) 
                self['QC' ] = self['QC' ].str.cat(second_slice[:, 1])

        self.__class__ = fastqDF
        return self.drop_abnormal_lengths() if enforce else self

    def expected_errors(self, maxEE=None):
        try:
            EE = self['QC'].apply(_expected_errors)
        except ValueError as e:
            print(self.short_filename, 'had an invalid Q score.')
            raise e
        if maxEE is None:
            return EE
        droped = self.loc[EE <= maxEE, :]
        droped.__class__ = fastqDF
        EE_total = EE.sum()
        return droped, EE_total

    def write(self, filename):
        if '__plus_sign__' not in self.columns:
            self.insert(2, '__plus_sign__', len(self)*['+'])
        filename = filename.split('.fastq')[0] + '.fastq'
        self.to_csv(filename + '.gz', 
                    columns=self.essential_columns[0:2] + ['__plus_sign__'] + self.essential_columns[2:], 
                    compression='gzip',
                    header=False,
                    index=False,
                    sep='\n')
            
class logPrint(object):
    def __init__(self, filename):
        import sys
        self.filename = filename 
        self.f = open(filename, 'w')
        self.f.write("# Output Summary of {:}:\n".format(*sys.argv[:1]))
        self.f.write("# --------------------------------------\n")

    def __call__(self, line):
        print(line)
        self.f.write('# '+line+'\n')        

    def close(self):
        self.f.close()

NW_kwargs = dict(gap_open=-10, gap_extend=-0.2, matrix='/home/chris/Dropbox/winslow/barcode_analysis/CM1')

ascii_gap = ord('-')
ascii_N = ord("N")

def repair_N(seq, ref):
    from nwalign import global_align
    ref_align, seq_align = global_align(ref, seq, **NW_kwargs)
    i = 0
    for r, s in zip(ref_align, seq_align):
        if s == ascii_N and r != ascii_gap:
            seq[i] = r
        i += s != ascii_gap
    return seq.tostring().decode("ascii")

def NW_fit(seq, ref, flanking_nucleotides_for_scoring):
    """score=score of [match, mismatch, gap_start, gap_extend]"""
    import nwalign
    ref_align, seq_align = nwalign.global_align(ref, seq, **NW_kwargs)
    ref_align.reverse()
    head = seq_align[:ref_align.index(ascii_N)]
    tail_start = len(ref_align) - ref_align.index(ascii_N)
    tail = seq_align[tail_start:]
    
    start = len(head) - head.count(ascii_gap)
    stop = len(seq_align) - len(tail) + tail.count(ascii_gap)
    score_start = max(start - flanking_nucleotides_for_scoring, 0)
    score_stop = stop + flanking_nucleotides_for_scoring
    seq2 = seq[score_start:score_stop]
    ref2 = ref[score_start:score_stop] 
    score = nwalign.score_alignment(ref2, seq2, **NW_kwargs)
    return start, stop, score

import re
class singleMismatcher(object):
    def __init__(self, substring):
        self.exact_match = substring
        self.patterns = [substring[:i]+'.'+substring[i+1:] for i in range(0, len(substring)-1)] + [substring[:-1]]
        self.re_objects = list(map(re.compile, self.patterns))

    def find(self, searchstring):
        for re_object in self.re_objects:
            searched = re_object.search(searchstring)
            if searched is not None:
                return searched.start()
        else:
            return -1

def hamming_distance(a, b):
    from scipy.spatial.distance import hamming
    N = len(a)
    return round(N*hamming(np.fromiter(a, 'S1', N), np.fromiter(b, 'S1', N)))

cdef double QC_map_C[256] 
cdef double [:] QC_map = QC_map_C
cdef double dp = 0.7943282347
QC_map[33] = 1.0
for i in range(34, 76):
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

