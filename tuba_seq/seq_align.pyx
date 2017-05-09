"""
A simple wrapper around Isaac Turner's seq-align repo:

Smith-Waterman & Needleman-Wunsch Alignment in C
url: https://github.com/noporpoise/seq-align
author: Isaac Turner turner.isaac@gmail.com
license: Public Domain

The python class NW performs Needleman-Wunsch Alignments with keyword-argument
handling of the scoring matrix & alignment parameters, and assimilates memory 
management--it does not do anything else. 
"""

cdef extern from 'alignment.h':
    ctypedef struct aligner_t:
        pass

ctypedef aligner_t nw_aligner_t

cdef extern from "needleman_wunsch.h":
    ctypedef struct alignment_t:
        char *result_a
        char *result_b
        int score

    void needleman_wunsch_free(nw_aligner_t *nw)

cdef extern from "needleman_wunsch.c":
    void needleman_wunsch_align(const char *a, const char *b,
                                const scoring_t *scoring,
                                nw_aligner_t *nw, alignment_t *result)
    nw_aligner_t* needleman_wunsch_new()

cdef extern from "alignment.c":
    alignment_t* alignment_create(size_t capacity)

cdef extern from "alignment_scoring.h":
    ctypedef struct scoring_t:
        pass

cdef extern from "alignment_scoring.c":
    void scoring_add_mutation(scoring_t* scoring, char a, char b, int score)
    
    void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  int no_start_gap_penalty, int no_end_gap_penalty,
                  int no_gaps_in_a, int no_gaps_in_b,
                  int no_mismatches, int case_sensitive)

cdef class NW(object):
    cdef nw_aligner_t *_aligner 
    cdef scoring_t _scoring

    def __cinit__(self):
        self._aligner = needleman_wunsch_new()  

    def __dealloc__(self):
        cdef nw_aligner_t *aligner = self._aligner
        needleman_wunsch_free(aligner)

    def __init__(self, int match=1, int mismatch=-2, int gap_open=-4, int gap_extend=-1,
                    bint start_gap_penalty=False, bint end_gap_penalty=False, 
                    bint gaps_in_a=True, bint gaps_in_b=True, bint mismatches=True, bint case_sensitive=True):
        
        # I do not understand why the boolean variables were negated.
        # Also, cython's bint type was exhibiting peculiar behavior, so I declared the bools explicitly. 

        cdef:
            int no_start_gap_penalty = 0 if start_gap_penalty else 1
            int no_end_gap_penalty = 0 if end_gap_penalty else 1
            int no_gaps_in_a = 0 if gaps_in_a else 1
            int no_gaps_in_b = 0 if gaps_in_b else 1
            int no_mismatches = 0 if mismatches else 1 
            int c_case_sensitive = 1 if case_sensitive else 0
        cdef scoring_t *scoring = &(self._scoring)

        scoring_init(scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, c_case_sensitive)
    
    cpdef add_mutation(self, a, b, int score):
        char_a = a.encode('ascii')
        char_b = b.encode('ascii')
        cdef scoring_t *scoring = &(self._scoring)
        scoring_add_mutation(scoring, char_a[0], char_b[0], score)

    cpdef add_neutral_N(self):
        cdef scoring_t *scoring = &(self._scoring)
        scoring_add_mutation(scoring, 'N', 'N', 0)
        
        scoring_add_mutation(scoring, 'N', 'A', 0)
        scoring_add_mutation(scoring, 'N', 'C', 0)
        scoring_add_mutation(scoring, 'N', 'G', 0)
        scoring_add_mutation(scoring, 'N', 'T', 0)
        
        scoring_add_mutation(scoring, 'A', 'N', 0)
        scoring_add_mutation(scoring, 'C', 'N', 0)
        scoring_add_mutation(scoring, 'G', 'N', 0)
        scoring_add_mutation(scoring, 'T', 'N', 0)

    cpdef char_align(self, char *seq_a, char *seq_b):
        cdef size_t len_a = len(seq_a)
        cdef size_t len_b = len(seq_b)
        cdef alignment_t *result = alignment_create(len_a + len_b)
        cdef scoring_t *scoring = &(self._scoring)
        cdef nw_aligner_t *aligner = self._aligner
        needleman_wunsch_align(seq_a, seq_b, scoring, aligner, result)
        return result.result_a, result.result_b, result.score
    
    cpdef align(self, a, b):
        r_a, r_b, score = self.char_align(a.encode('ascii'), b.encode('ascii'))
        return r_a.decode('ascii'), r_b.decode('ascii'), score
   
    cpdef char_score(self, char *seq_a, char *seq_b):
        cdef size_t len_a = len(seq_a)
        cdef size_t len_b = len(seq_b)
        cdef alignment_t *result = alignment_create(len_a + len_b)
        cdef scoring_t *scoring = &(self._scoring)
        cdef nw_aligner_t *aligner = self._aligner
        needleman_wunsch_align(seq_a, seq_b, scoring, aligner, result)
        return result.score
    
