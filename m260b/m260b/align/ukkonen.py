"""\
Library for Ukkonen algorithm based local alignment

"""
import os

import numpy as np

from m260b.align.c_smithwaterman import banded_sw_c

# smith watrman params
MATCH_GAIN, MATCH_GAIN_LENIENT = 15, 30
MISMATCH_COST, MISMATCH_COST_LENIENT = 9, 12
DEL_GAP_OPEN_COST, DEL_GAP_OPEN_COST_LENIENT = 20, 10
INS_GAP_OPEN_COST, INS_GAP_OPEN_COST_LENIENT = 9, 5
GAP_EXTENSION_COST, GAP_EXTENSION_COST_LENIENT = 3, 0
MOV_TO_CIG = {0: 'M', 1: 'D', -1: 'I'}

BAND_DIFF = 100
BAND_START_BEST = 100

def banded_sw(ref, read, not_in_ref_penalty=10):
    gain, moves = _banded_sw_matrix(ref, read)
    return produce_alignment(ref, read, gain, moves)


def full_sw(ref, read, ref_cost=10, lenient=False):
    gain, moves = _full_sw_matrix(ref, read, ref_cost, lenient)
    return produce_alignment(ref, read, gain, moves)


_SW_SCORE, _SW_MOVE = None, None
def _banded_sw_matrix(ref, read, not_in_ref_penalty=10):
    global _SW_SCORE, _SW_MOVE
    if _SW_SCORE is None or _SW_SCORE.shape != (len(ref) + 1, len(read) + 1):
        _SW_SCORE = np.zeros((len(ref) + 1, len(read) + 1), dtype=np.int32)
        _SW_MOVE = np.zeros((len(ref) + 1, len(read) + 1), dtype=np.int32)
    # score and move are now cached
    return banded_sw_c(str(ref), str(read), len(ref), len(read), not_in_ref_penalty, _SW_SCORE, _SW_MOVE) 


def _full_sw_matrix(ref, read, ref_cost=10, lenient=False):
    """\
    Implementation of the full smith-waterman algorithm. A prelude to banded smith-waterman
    and eventually Ukkonen.

    """
    swgain = np.zeros((1+len(ref), 1+len(read)), dtype=int)  # holds the SW gain value
    moves = np.zeros((1+len(ref), 1+len(read)), dtype=int)   # holds the move for each position
    swgain[0, 1:] = -ref_cost
    mgain, mcost, dgap, igap, gapext = MATCH_GAIN, MISMATCH_COST, DEL_GAP_OPEN_COST, INS_GAP_OPEN_COST, \
                                         GAP_EXTENSION_COST
    if lenient:
        mgain, mcost, dgap, igap, gapext = MATCH_GAIN_LENIENT, MISMATCH_COST_LENIENT, \
                                             DEL_GAP_OPEN_COST_LENIENT, INS_GAP_OPEN_COST_LENIENT, \
                                             GAP_EXTENSION_COST_LENIENT
    # 0 = diagonal, +1 = vertical, -1 = horizontal
    # fill in the matrix
    for row in (1 + r for r in xrange(len(ref))):
        best_row = - np.inf
        for col in (1 + c for c in xrange(len(read))):
            if ref[row-1] == read[col-1]:
                Mvalue = swgain[row-1, col-1] + mgain 
            else:
                Mvalue = swgain[row-1, col-1] - mcost 
            if moves[row-1, col] == 1:
                Dvalue = swgain[row-1, col] - gapext 
            else:
                Dvalue = swgain[row-1, col] - dgap 
            if moves[row, col-1] == -1:
                Ivalue = swgain[row, col-1] - gapext 
            else:
                Ivalue = swgain[row, col-1] - igap 
            best = max(Mvalue, Dvalue, Ivalue)
            swgain[row, col] = best
            if Ivalue == best:
                moves[row, col] = -1
            elif Dvalue == best:
                moves[row, col] = 1
            if best > best_row:
                best_row = best
    return swgain, moves 


def produce_alignment(ref, read, swgain, moves):
    # since the read is shorter than the ref, get the best row for the last column
    if swgain is None:
        return None, None, None, None
    row, col = np.argmax(swgain[:, -1]), len(read)
    score = swgain[row, col]
    best_moves, last_move, n_same_move, mismatches = [], None, 0, 0
    while row > 0 and col > 0:  # remember row 0 is just the 'string start position' row
         if moves[row, col] == last_move:
             n_same_move += 1
         else:
             if last_move is not None:
                 best_moves.append((last_move, n_same_move))
             last_move = moves[row, col]
             n_same_move = 1
         if moves[row, col] == 0:
             if ref[row-1] != read[col-1]:
                 mismatches += 1
             row, col = row - 1, col - 1
         elif moves[row, col] == 1:
             row = row - 1
         else:
             col = col - 1
    best_moves.append((last_move, n_same_move))
    # get the relative position of the read to the reference
    offset_into_ref = row - col
    # convert the moves directly into a cigar
    cigar = None
    try:
        cigar = [[MOV_TO_CIG[i], j] for i, j in best_moves][::-1]
    except KeyError:
        raise ValueError('Error in smith-waterman\nref={}\nread={}\nbest_moves={}\nscores=\n{}'.format(ref, read, best_moves, swgain))
    # cannot start with an insertion
    if cigar[0][0] == 'I':
        cigar[0][0] = 'S'
    # cannot end either
    if cigar[-1][0] == 'I':
        cigar[-1][0] = 'S'
    return offset_into_ref, ''.join(['{}{}'.format(e[1], e[0]) for e in cigar]), score, mismatches 

"""
def explore_row(row, costbound, lss_matrix, max_cost):
    ""
    Fill out a row of the longest substring matrix. This is done with the following
    boundary conditions:
               
     (boundary condition)       
       row    V  0  1  2  3  4  5   cost
         0  [-1][ ][ ][ ][ ][ ][ ]
         1      [K][ ][ ][ ][ ][ ]
         2      [K][K][o][ ][ ][ ]
         3      [K][K][K][ ][ ][ ]
         4      [K][K][K][K][ ][ ]
    ""



def ukkonen_aln(ref, read, max_cost=6):
    ""\
    Use Ukknonen's algorithm to locally align the read sequence to the reference
    sequence. This works just with the bases and is *not* quality score aware.

    ""
    ref_pos, read_pos, distance, cost_bound = 0, 0, 0, 0
    last_diag = len(ref)-len(read)  # the last of the diagonal in the str x str matrix
    main_diagonal = (max_cost - last_diag) / 2  # the main diagonal of the SW matrix shows up as this row
    lss_mat = - np.ones(diag, max_cost) # initialization, set longest string size to -1 
    while distance != len(ref):
         # get the insert cost
         ins_len = explore_row(last_diag + 1, cost_bound - GAP_COST, lss_mat, max_cost)
         # get the deletion cost
         del_len = explore_row(last_diag - 1, cost_bound - GAP_COST, lss_mat, max_cost) + 1
         # and the mismatch cost
         mis_len = explore_row(last_diag, cost_bound - MISMATCH_COST, lss_mat, max_cost) + 1
         # the cost of this entry is the largest of the above
         longest_string = max(ins_len, del_len, mis_len)
         # note that since the row is a diagonal of the SW matrix, we can index
         # into the strings using the row and longest_string as coordinates
         if longest_string < 0:
             return None  # could not align
         while ref[longest_string + 1] == read[longest_string - last_diag + 1]:
             longest_string += 1
         lss_mat[last_diag, cost_bound] = longest_string
         distance = longest_string
         cost_bound += 1
    # get the longest matching prefix
    cost_mat[0, 0] = len(os.path.commonprefix([ref, read]))  # such utilities! much good!
    if diag_len > distance:
        cost_mat[diag_len, distance] = - np.inf
    while cost < max_cost:
"""    
