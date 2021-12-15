#!/usr/bin/env python3
"""
Simple python wrapper for SSW library
Please put the path of libssw.so into LD_LIBRARY_PATH or pass it explicitly as a parameter
By Yongan Zhao (March 2016)
"""

import ctypes as ct
import numpy as np
from tuba_seq.ssw_lib import CSsw

from pathlib import Path
ssw_lib = Path.home() / 'tuba-seq/Complete-Striped-Smith-Waterman-Library/src'
ssw = CSsw(str(ssw_lib))

def construct_alignment(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = b'MIDNSHP=X'
    sCigar = b''
    sQ = []
    sR = []
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        c = 77 if m > 8 else sCigarInfo[m]
        if c == 77:  #'M'
            sQ.append(q[nQOff : nQOff+n])
            sR.append(r[nROff : nROff+n])
            nQOff += n
            nROff += n
        elif c == 73: #'I'
            sQ.append(q[nQOff : nQOff+n])
            sR.append(b'-' * n)
            nQOff += n
        elif c == 68: #'D'
            sQ.append(b'-' * n)
            sR.append(r[nROff : nROff+n])
            nROff += n
        else:
            raise ValueError("Invalid Cigar Annotation ({:})".format(c))
    return b''.join(sQ), b''.join(sR)

o_N = ord(b'N')
o_gap = ord(b'-')

class SW(object):
    def __dealloc__(self):
        ssw.init_destroy(self.qProfile)
        ssw.align_destroy(self.res)

    def __init__(self, match=6, mismatch=2, gap_open=6, gap_extend=1):
        # init DNA score matrix
        lEle = list(b'ACGTN')
        nInt2Ele = np.array(lEle, dtype='int8')
        self.nEle2Int = np.zeros(256, dtype='int8')
        self.nEle2Int[nInt2Ele] = np.arange(len(lEle))
        self.lScore = [0 if nuc1 == o_N or nuc2 == o_N else (match if nuc1 == nuc2 else mismatch) for nuc1 in lEle for nuc2 in lEle]
            # translate score matrix to ctypes
        self.mat = (len(self.lScore) * ct.c_int8) ()
        self.mat[:] = self.lScore
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.lEle = lEle

    def char_align(self, char_query, char_reference):
        nQuery = self.nEle2Int[list(char_query)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.qProfile = ssw.ssw_init(nQuery, ct.c_int32(len(char_query)), self.mat, len(self.lEle), 2)
        nMaskLen = len(char_query) // 2 if len(char_query) > 30 else 15
        nFlag = 2
        nReference = self.nEle2Int[list(char_reference)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.res = ssw.ssw_align(self.qProfile, nReference, ct.c_int32(len(char_reference)), self.gap_open, self.gap_extend, nFlag, 0, 0, nMaskLen)
        contents = self.res.contents
        lCigar = [contents.sCigar[idx] for idx in range(contents.nCigarLen)]
        out = construct_alignment(char_query, char_reference, contents.nQryBeg, contents.nRefBeg, lCigar)
        return out

    def find_N_start_stop(self, char_query, char_reference):
        query_align, ref_align = self.char_align(char_query, char_reference)
        query_begin = self.res.contents.nQryBeg
        N_start = ref_align.find(b'N')
        N_stop = ref_align.rfind(b'N') + 1
        head_gaps = query_align[0:N_start].count(b'-')
        barcode_gaps = query_align[N_start:N_stop].count(b'-')
        return query_begin + N_start - head_gaps, query_begin + N_stop - head_gaps - barcode_gaps

    def ClonTracer_start(self, char_query, char_reference):
        query_align, ref_align = self.char_align(char_query, char_reference)
        return self.res.contents.nQryBeg

    def fill_Ns(self, char_query, char_reference):
        query_align, ref_align = self.char_align(char_query, char_reference)
        middle = bytes(bytearray([q if (q != o_N or r == o_gap) else r for q, r in zip(query_align, ref_align) if q != o_gap]))
        contents = self.res.contents
        out = char_query[:contents.nQryBeg] + middle + char_query[contents.nQryEnd+1:]
        return out

    def align(self, query, reference):
        return self.char_align(query.encode('ascii'), reference.encode('ascii'))

    def char_score(self, char_query, char_reference):
        nQuery = self.nEle2Int[list(char_query)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.qProfile = ssw.ssw_init(nQuery, ct.c_int32(len(char_query)), self.mat, len(self.lEle), 2)
        nMaskLen = len(char_query) // 2 if len(char_query) > 30 else 15
        nFlag = 0
        nReference = self.nEle2Int[list(char_reference)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.res = ssw.ssw_align(self.qProfile, nReference, ct.c_int32(len(char_reference)), self.gap_open, self.gap_extend, nFlag, 0, 0, nMaskLen)
        return self.res.contents.nScore

    def fill_Ns(self, char_query, char_reference):
        query_align, ref_align = self.char_align(char_query, char_reference)
        middle = bytes(bytearray([q if (q != o_N or r == o_gap) else r for q, r in zip(query_align, ref_align) if q != o_gap]))
        contents = self.res.contents
        out = char_query[:contents.nQryBeg] + middle + char_query[contents.nQryEnd+1:]
        return out

    def align(self, query, reference):
        return self.char_align(query.encode('ascii'), reference.encode('ascii'))

    def char_score(self, char_query, char_reference):
        nQuery = self.nEle2Int[list(char_query)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.qProfile = ssw.ssw_init(nQuery, ct.c_int32(len(char_query)), self.mat, len(self.lEle), 2)
        nMaskLen = len(char_query) // 2 if len(char_query) > 30 else 15
        nFlag = 0
        nReference = self.nEle2Int[list(char_reference)].ctypes.data_as(ct.POINTER(ct.c_int8))
        self.res = ssw.ssw_align(self.qProfile, nReference, ct.c_int32(len(char_reference)), self.gap_open, self.gap_extend, nFlag, 0, 0, nMaskLen)
        return self.res.contents.nScore
