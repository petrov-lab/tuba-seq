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
        data = pd.read_table(filename, names=['__'], encoding='ascii', quoting=3, dtype=str)['__'].values.reshape((-1, 4))[:, [True, True, False, True]] 
        self = cls(pd.DataFrame(data=data,
                    columns=cls.essential_columns))
        if use_Illumina_filter:
            self.select_reads(self['header'].str.contains(":N:"))       # Y = filtered, N = passed filter
       
        sample_header = self.loc[0, 'header']
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
c_gap = '-'.encode('ascii') 
c_N = 'N'.encode('ascii') 

def repair_N(c_seq, c_ref):
    seq_align, ref_align, _score = nw.char_align(c_seq, c_ref)
    seq_align, ref_align = seq_align.decode("ascii"), ref_align.decode("ascii")
    return ''.join([r_i if (s_i == 'N' and r_i != '-') else s_i for i, (s_i, r_i) in enumerate(zip(seq_align, ref_align)) if s_i != '-'])

def find_start_stop(seq, ref): 
    """score=score of [match, mismatch, gap_start, gap_extend]"""
    seq_align, ref_align, _score = nw.char_align(seq, ref) 
    cdef:
        int ref_start = ref_align.index(c_N)
        int ref_stop = ref_align.rindex(c_N) + 1
        int head_gaps = seq_align[0:ref_start].count(c_gap)
        int barcode_gaps = seq_align[ref_start:ref_stop].count(c_gap)
    return ref_start - head_gaps, ref_stop - head_gaps - barcode_gaps

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

