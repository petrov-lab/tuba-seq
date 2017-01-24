import numpy as np
cimport numpy as np

import array
from cpython cimport array

cimport cython
import sys
import os.path

cdef extern from "stdlib.h":
    ctypedef unsigned int size_t
    size_t strlen(char *s)


ctypedef np.int_t DTYPE_INT
ctypedef np.uint_t DTYPE_UINT
ctypedef np.int8_t DTYPE_BOOL

cdef size_t UP = 1, LEFT = 2, DIAG = 3, NONE = 4

cdef inline int imax2(int a, int b):
    if a >= b: return a
    return b

#def score_alignment(unsigned chara, b, int gap_open, int gap_extend, matrix):
def score_alignment(unsigned char[:] al, unsigned char[:] bl, int gap_open, int gap_extend, matrix):
#    cdef char *al = a
#    cdef char *bl = b
    cdef size_t l = al.shape[0], i 
    cdef int score = 0, this_score
    assert bl.shape[0] == l, "alignment lengths must be the same"
    cdef np.ndarray[DTYPE_INT, ndim=2] mat
    mat = read_matrix(matrix)

    cdef bint gap_started = 0

    for i in range(l):
        if al[i] == c"-" or bl[i] == c"-":
            score += gap_extend if gap_started else gap_open
            gap_started = 1
        else:
            this_score = mat[al[i], bl[i]]
            score += this_score
            gap_started = 0
    return score

cdef read_matrix(path, dict cache={}):
    """
    so here, we read a matrix in the NCBI format and put
    it into a numpy array. so the score for a 'C' changing
    to an 'A' is stored in the matrix as:
        mat[ord('C'), ord('A')] = score
    as such, it's a direct array lookup from each pair in the alignment
    to a score. this makes if very fast. the cost is in terms of space.
    though it's usually less than 100*100.
    """
    if path in cache: 
        return cache[path]
    cdef np.ndarray[DTYPE_INT, ndim=2] a
    cdef size_t ai = 0, i
    cdef int v, mat_size

    fh = open(path)
    headers = None
    while headers is None:
        line = fh.readline().strip()
        if line[0] == '#': continue
        headers = [ord(x) for x in line.split(' ') if x]
    mat_size = max(headers) + 1

    a = np.zeros((mat_size, mat_size), dtype=np.int)

    line = fh.readline()
    while line:
        line_vals = [int(x) for x in line[:-1].split(' ')[1:] if x]
        for ohidx, val in zip(headers, line_vals):
            a[headers[ai], ohidx] = val
        ai += 1
        line = fh.readline()
        
    cache[path] = a
    return a

#def global_align(object _seqj, object _seqi, int gap_open=-1, int gap_extend=-1, object matrix):
@cython.boundscheck(False)
@cython.nonecheck(False)
def global_align(unsigned char[:] seqj, unsigned char[:] seqi, int gap_open=-1, int gap_extend=-1, object matrix='CM1'):
    """
    perform a global sequence alignment (needleman-wunsch) on seq and and 2. using
    the matrix for nucleotide transition from matrix if available.
    where matrix is of the format provided in the ncbi/data directory.

    >>> from nwalign import global_align
    >>> global_align('COELANCANTH', 'PELICAN')
    ('COELANCANTH', '-PEL-ICAN--')

    """

    cdef bint flip = 0

#    cdef char* seqj = _seqj
#    cdef char* seqi = _seqi

    cdef size_t max_j = seqj.shape[0]
    cdef size_t max_i = seqi.shape[0] 
    if max_i == max_j == 0:
        return "", ""

    cdef np.ndarray[DTYPE_INT, ndim=2] amatrix
    amatrix = read_matrix(matrix)

    if max_j > max_i:
        flip = 1
        seqi, seqj = seqj, seqi
        max_i, max_j = max_j, max_i

    # need to initialize j for the case when it's a zero-length string.
    #cdef size_t i = 0, j = 0, seqlen, align_counter = 0, p
    cdef size_t i = 0, j = 0, seqlen, p
    cdef int diag_score, up_score, left_score, tscore

    #cdef char *align_j, *align_i
    cdef char ci, cj
    #cdef PyObject *ai, *aj
    cdef int zero=0, one=1

    assert gap_extend <= 0, "gap_extend penalty must be <= 0"
    assert gap_open <= 0, "gap_open must be <= 0"

    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_i = np.ones((max_i + 1,), dtype=np.int8)
    cdef np.ndarray[DTYPE_BOOL, ndim=1] agap_j = np.ones((max_j + 1,), dtype=np.int8)

    cdef np.ndarray[DTYPE_INT, ndim=2] score = np.empty((max_i + 1, max_j + 1), dtype=np.int)
    cdef np.ndarray[DTYPE_UINT, ndim=2] pointer = np.empty((max_i + 1, max_j + 1), dtype=np.uint)
  
    pointer[0, 0] = NONE
    score[0, 0] = 0
    
    pointer[0, 1:] = LEFT
    pointer[1:, 0] = UP

    score[0, 1:] = gap_open + gap_extend * np.arange(0, max_j, dtype=np.int)
    score[1:, 0] = gap_open + gap_extend * np.arange(0, max_i, dtype=np.int)

    agap_i[0] = zero

    for i in range(1, max_i + 1):
        ci = seqi[<size_t>(i - 1)]

        agap_j[0] = zero
        for j in range(1, max_j + 1):
            agap_j[j] = one
            cj = seqj[<size_t>(j - 1)]
            tscore = amatrix[ci, cj]

            diag_score = score[<size_t>(i - 1), <size_t>(j - 1)] + tscore
            up_score   = score[<size_t>(i - 1), j] + (gap_open if agap_i[<size_t>(i - 1)] == zero else gap_extend)
            left_score = score[i, <size_t>(j - 1)] + (gap_open if agap_j[<size_t>(j - 1)] == zero else gap_extend)

            """
            fix cases where scores are tied.
            choose diagonal when not at the ends of either string.
            and choose up/left (gap) when at the end of the string.
            this forces gaps to the ends of the aligments.
            """
            if diag_score == left_score:
                # so here. if we're at the end we choose a gap.
                # otherwise, we choose diag
                if i == max_i or i == 1:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            elif diag_score == up_score:
                if j == max_j or j == 1:
                    score[i, j] = up_score
                    pointer[i, j] = UP
                else: # we want a diagonal
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero
            # end of ambiguous score checks.

            elif up_score > diag_score: #checked.
                if up_score > left_score:
                    score[i, j]  = up_score
                    pointer[i, j] = UP
                else:
                    score[i, j]   = left_score
                    pointer[i, j] = LEFT
            elif diag_score > up_score:
                if left_score > diag_score:
                    score[i, j] = left_score
                    pointer[i, j] = LEFT
                else:
                    score[i, j] = diag_score
                    pointer[i, j] = DIAG
                    agap_i[i] = zero

    seqlen = max_i + max_j
    #ai = PyString_FromStringAndSize(NULL, seqlen)
    #aj = PyString_FromStringAndSize(NULL, seqlen)

    cdef int score_max #, = score[:, -1].max()

    # had to use this and PyObject instead of assigning directly...
    #align_j = PyString_AS_STRING(aj)
    #align_i = PyString_AS_STRING(ai)

    cdef array.array char_array_template = array.array('B', []) #, gap_str
    cdef array.array align_i = array.clone(char_array_template, seqlen, False)
    cdef array.array align_j = array.clone(char_array_template, seqlen, False)
    cdef size_t align_counter = seqlen - 1

    p = pointer[i, j]
    while p != NONE:
        if p == DIAG:
            i -= 1
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = seqi[i]
        elif p == LEFT:
            j -= 1
            align_j[align_counter] = seqj[j]
            align_i[align_counter] = c"-"
        elif p == UP:
            i -= 1
            align_j[align_counter] = c"-"
            align_i[align_counter] = seqi[i]
        else:
            raise Exception('wtf!:pointer: %i', p)
        #align_counter += 1
        align_counter -= 1
        p = pointer[i, j]


    #_PyString_Resize(&aj, align_counter)
    #_PyString_Resize(&ai, align_counter)
    #if flip:
    #    return (<object>ai)[::-1], (<object>aj)[::-1] #, score.max()
    #else:
    #    return (<object>aj)[::-1], (<object>ai)[::-1] #, score.max()
    align_counter += 1
    #return (align_i[align_counter:].tostring().decode('ascii'), align_j[align_counter:].tostring().decode('ascii') ) if flip else (align_j[align_counter:].tostring().decode('ascii'), align_i[align_counter:].tostring().decode('ascii') )
    return (align_i[align_counter:].tostring(), align_j[align_counter:]) if flip else (align_j[align_counter:], align_i[align_counter:] )
    #return (align_i[align_counter:].tostring(), align_j[align_counter:].tostring() ) if flip else (align_j[align_counter:].tostring(), align_i[align_counter:].tostring() )

