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

c_gap = '-'.encode('ascii') 
c_N = 'N'.encode('ascii') 

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
        
        self.barcode_length = len(self.ref) - 2*self.alignment_flank
        
        self.c_ref = self.ref.encode('ascii')
        self.max_score = nw.char_score(self.c_ref, self.c_ref)
        self.min_align_score = args.min_align_score

        aft_length = len(self.full) - stop
        pre_slice = (slice(start-aft_length, stop+aft_length) if start > aft_length  else slice(args.trim, stop+start-args.trim)) if args.symmetric_flanks else slice(0, len(self.full)) 
        self.pre_slice = slice(pre_slice.start + args.trim, pre_slice.stop - args.trim)
        
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
        return ''.join([s if (s != 'N' or r == '-') else r for s, r in zip(seq_align.decode('ascii'), ref_align.decode('ascii')) if s != '-'])
    
    def iter_fastq(self, sample, filenames, input_file):
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

        LINE_3 = '\n+\n'.encode('ascii')
        END = '\n'.encode('ascii')
        with open(filenames[0], 'wb') as training_file, open(filenames[1], 'wb') as cluster_file:
            while True:
                header = input_file.readline()
                if not header:
                    break
                c_DNA = input_file.readline()[self.pre_slice]
                input_file.readline()
                QC = input_file.readline()[self.pre_slice]
                if not QC:
                    raise RuntimeError("Input FASTQ file was not 4x lines long")
                DNA = c_DNA.decode('ascii')
                if ':Y:'.encode('ascii') in header:
                    Filtered += 1
                    continue
                if DNA in self.alignments:
                    start, stop, score = self.alignments[DNA]
                else:
                    start, stop = self.find_start_stop(c_DNA)
                    c_DNA_scoring = c_DNA[start - self.alignment_flank:stop + self.alignment_flank]
                    score = nw.char_score(c_DNA_scoring, self.c_ref)/self.max_score
                    self.alignments[DNA] = start, stop, score
                scores[score] += 1
                if score < self.min_align_score:
                    Unaligned += 1
                    self.unaligned[DNA] = 1 + self.unaligned.get(DNA, 0)
                    continue
                if abs((stop - start) - BL) > self.allowable_deviation:
                    Wrong_Barcode_Length += 1
                    continue
                training_DNA = DNA[start - TF:start]+DNA[stop:stop + TF]
                if 'N' not in training_DNA and len(training_DNA) == 2*TF:
                    tQC = QC[start-TF:start]+QC[stop:stop+TF]
                    assert len(tQC) == len(training_DNA), '{:}\n{:}'.format(training_DNA, tQC)
                    training_file.write(header+training_DNA.encode('ascii')+LINE_3+tQC+END)
                if 'N' in DNA:
                    DNA = self.repair_N(c_DNA)
                cluster_DNA = DNA[start - CF:start+BL+CF]
                if 'N' in cluster_DNA:
                    Residual_N += 1
                    continue
                if len(cluster_DNA) != (BL+ 2*CF):
                    Insufficient_Flank += 1
                else:
                    Clustered += 1
                    cQC = QC[start-CF:start+BL+CF]
                    assert len(cQC) == len(cluster_DNA), '{:}\n{:}'.format(cluster_DNA, cQC)
                    cluster_file.write(header+cluster_DNA.encode('ascii')+LINE_3+cQC+END)

        self.instruments[sample] = header.decode('ascii').split(':')[0]
        self.scores[sample] = pd.Series(scores)
        return pd.Series([Filtered,   Unaligned,   Wrong_Barcode_Length,   Residual_N,   Insufficient_Flank,   Clustered], 
         index=pd.Index(['Filtered', 'Unaligned', 'Wrong Barcode Length', 'Residual N', 'Insufficient Flank', 'Clustered']))

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

