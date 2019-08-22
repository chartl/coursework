"""\
Instead of traversing locus by locus and forming pileups, instead collect reads into windows,
check to see if there is something going on there, and if so, return the read pool (or
appropriate pileups, if desired).

"""
from collections import namedtuple
from itertools import islice, izip

import numpy as np

from m260b.traversal.pileup import make_pileups


CIGAR_TO_SCORE = {0: 0, # match/mismatch
                  1: 4, # insertion
                  2: 3, # deletion
                  4: 1} # soft clip
Region = namedtuple('Region', ('start', 'stop', 'reference'))


def sliding_window(lst, size):
   """\
   Generate a sliding window of size `size` over the list.

   This is straight out of itertools documentation

   """
   it = iter(lst)
   result = list(islice(it, size))
   if len(result) == size:
       yield result
   for elem in it:
       result = result[1:] + [elem]
       yield result


def active_pileups(reads, ref, chr, start, flank=10, dfrac=0.5, pwin=0):
    """\
    Use the active region traversal to generate pileups

    """
    if pwin > 2*flank:
        raise ValueError('Need flank of size at least {} for pwin {}'.format(pwin/2, pwin))
    regions = active_regions(reads, ref, chr, start, flank, dfrac)
    for region, reads in regions:
        for pileup in make_pileups(reads, ref, chr, offset=region.start, window_size=pwin, stop=region.stop):
            yield pileup 


def _end_offset(read):
    if read.cigar[-1][0] != 4:
        return read.aend - 1
    return read.aend - 1 + read.cigar[-1][1]  # aend disregards soft clip. I don't want to.



def active_regions(reads, reference, contig, start_offset, flank=10, dfrac=0.5):
    """\
    Collects reads into a pool, checks if the pool is active, extends it if need be, and
    returns an iterator over the reads in the pool.

    The efficiency gain here is that most regions will be inactive, and in order to
    check that a region is active, we only need to use some fraction (`dfrac`) of the
    reads.

    Let's say we want to check a certain location in the reference
                          v
        ------------------^--------------------

    Then we open up a window of `flank` bps and collect reads
                          v
        -------|----------^----------|-----------
                       <------
            ------->  ------->   <------
               ------>     ------>
                   <------          ------->

    Now we `mark` the reference with activity scores
                          v      x
        -------|001000100082100094100|----------

    If it looks like there's another `active` site in the window
    open up another `flank`bp
                          v      x
        -------|----------^------^-------------|
           [reads reads reads reads reads reads reads]

    and repeat until there's not an active mark within `flank` of the
    window.

    If any position in the reference is active, return a generator over
    the read pool. Otherwise, continue looking.

    There are lots of intervals here, so care needs to be taken. The purpose
    of the flank is to include that many basepairs. So if we want the center
    of our window at 15 with flank size 8, then the interval is 
       [7, 23] INCLUSIVE for a total size of 23-7+1 = 17
    
    """
    rgen = np.random.RandomState(20141022)
    reads = iter(reads)  # in case not already an iterator
    refsize = len(reference)
    # advance the starting position to the first read overlapping
    # the start offset
    read = next(reads)
    while _end_offset(read) < start_offset:
        read = next(reads)
    # open up the window of 1 + 2*flank
    offset, midpt, win_stop= start_offset, start_offset + flank, start_offset + 2*flank
    pool = list()
    # create a reference score object
    refscore = np.zeros((win_stop-offset+1,), dtype=int)
    depth = np.zeros((win_stop-offset+1,), dtype=int)
    ref = reference[offset:win_stop]
    while win_stop < refsize: 
        # grab the reference seq
        ref = reference[offset:win_stop]
        while read.pos <= win_stop:  # reads that start on the right edge are fair game
            read_offset = read.pos - offset
            if read_offset > len(ref):
                raise ValueError('Error computing read offset from read\n{}\n'
                                 'offset: {}'.format(str(read), offset))
            # possibly update activity information (depending on fraction)
            if rgen.binomial(1, dfrac) == 1 and not read.is_unmapped:
                active_score = get_active_score(reference, read, offset, win_stop)
                _update_scores(refscore, depth, active_score, read_offset)
            pool.append(read)
            read = next(reads)
        is_active, extension_size = get_activity(refscore, depth, flank)
        if extension_size > 0:
            win_stop += extension_size  # will accumulate more reads
            refscore = np.hstack((refscore, np.zeros((extension_size,), dtype=int)))
            depth = np.hstack((depth, np.zeros((extension_size,), dtype=int)))
        else:
            if is_active:
                yield Region(offset, win_stop, ref), pool
            if win_stop == refsize - 1: 
                # done
                break
            else:
                pool, offset, win_stop = _advance_window(pool, offset, win_stop, flank, refsize)
                # TODO: do we need to remember which reads we used?
                refscore, depth = _advance_scores(pool, offset, win_stop, flank, 
                                                  reference, refscore, depth)


def _advance_scores(pool, offset, win_stop, flank, reference, refscore, depth):
    """\
    Shifting the window forward to its midpoint means we can reuse lots of the
    scores.
                  7               23
                  |-------o-------|
        score =  [abcdefghijklmnopq]


    Now we push forward
                 15
                  |-------o-------|
        score =  [ijklmnopq00000000]

    """
    #print(' '*35 + '{}|'.format(offset) + '-' * (flank-1) + 'o' + '-' * (flank-1) + '|{}'.format(win_stop))
    if refscore.shape[0] > (win_stop - offset + 1):
        refscore = refscore[-(win_stop-offset+1):]
        depth = depth[-(win_stop-offset+1):]
    new_start = refscore[flank:]
    refscore[:new_start.shape[0]] = new_start
    refscore[new_start.shape[0]:] = 0
    new_dp = depth[flank:]
    depth[:new_dp.shape[0]] = new_dp
    depth[new_dp.shape[0]:] = 0
    midpt = win_stop - flank
    # now push forward with the reads that we have
    for read in pool:
        # only care about things after the midpoint
        if _end_offset(read) > midpt:
            read_offset = read.pos - offset
            #print(' '*len(str(offset)) + ' '*(35 + read_offset) + '-'*(read.alen-1) + '>' + '   ({},{},{},{})'.format(read.pos, _end_offset(read), read.aend, read.cigarstring))
            _start_read = max(0, midpt - read.pos + 1) 
            _start_ref = midpt - offset + 1 
            active_score = get_active_score(reference, read, offset, win_stop, score_start=_start_read)
            refscore[_start_ref:(_start_ref + active_score.shape[0])] += active_score
            depth[_start_ref:(_start_ref + active_score.shape[0])] += 1
    return refscore, depth


def _update_scores(refscore, depth, active_score, read_offset):
    try:
        refscore[read_offset:(read_offset + active_score.shape[0])] += active_score
    except ValueError:
        raise ValueError('Error in attempting to update scores.\n'
                         'ref: {} (size={})\n'
                         'active: {}\n'
                         'offset: {}'.format(refscore.tolist(), len(refscore.tolist()),
                                             active_score.tolist(), read_offset))
    depth[read_offset:(read_offset + active_score.shape[0])] += 1


def _advance_window(pool, offset, win_stop, flank, refsize):
    newpool = list()
    for read in pool:
        if _end_offset(read) >= win_stop - flank:
            newpool.append(read)
        else:
            try:
                del _ACTIVE_SCORE_CACHE[_score_id(read)]
            except Exception:
                pass
    offset = win_stop - flank
    win_stop = min(offset + 2*flank, refsize - 1) 
    return newpool, offset, win_stop

def _score_id(read):
    return '{}.{}'.format(read.qname, 1 if read.is_read1 else 2)


_ACTIVE_SCORE_CACHE = dict()
def get_active_score(ref, read, win_start, win_end, score_start=None):
    """\
    Determine the activity score of the read from its mismatches and cigar elements.

    This will determine the full score of the read, and subset appropriately to
    the window size.

    """
    sstart = score_start or max(0, win_start - read.pos)
    send = win_end - _end_offset(read) - 1 if _end_offset(read) > win_end else read.alen 
    if _score_id(read) in _ACTIVE_SCORE_CACHE:
        return _ACTIVE_SCORE_CACHE[_score_id(read)][sstart:(1+send)]
    refpos, seqpos, score, carry = read.pos, 0, [], 0
    for idx, (oper, size) in enumerate(read.cigar):
        if refpos > _end_offset(read):
            raise ValueError('Should not be able to exceed read length')
        if oper == 0:  # match/mismatch
            next_ = [0 if ref[refpos + i] == b else 1
                     for i, b in enumerate(read.seq[seqpos:(seqpos+size)])]
            if carry:
                next_[0] += carry
                carry = 0 
            score.extend(next_)
            refpos += size
            seqpos += size
        elif oper == 1:  # insertion
            score[-1] += CIGAR_TO_SCORE[oper]
            seqpos += size
        elif oper == 2:  # deletion
            score[-1] += CIGAR_TO_SCORE[oper]
            score.extend([0 for _ in xrange(size)])
            if size > 1:
                score[-1] += CIGAR_TO_SCORE[oper]
            refpos += size
        elif oper == 4: # soft clip
            if idx is 0:
                # need a carryover
                carry = CIGAR_TO_SCORE[oper] * size
            else:
                score.extend([CIGAR_TO_SCORE[oper]] * size)
                refpos += size 
                seqpos += size
        else:
            raise ValueError('Unsupported cigar operator in {}'.format(read.cigarstring))
    score = np.array(score)
    _ACTIVE_SCORE_CACHE[_score_id(read)] = score
    if read.cigar[-1][0] != 4:
        assert score.shape[0] == read.alen, 'shape={}, cigar={}'.format(score.shape, read.cigarstring) 
    return score[sstart:(1+send)]


def get_activity(scores, depth, flank_size):
   """\
   Determine whether the activity scores for the given depth correspond to
   something fishy potentially happening in the sequence.

   Also if there's something fishy near the edge of the window, return the
   amount by which to extend the window.

   A window is active if any `3` positions have scores >= `20%` of their
   depth.

   """
   active = _get_active(scores, depth)
   is_active, extension = any(active), 0
   if is_active:
       activity = [i for i, t in enumerate(active) if t]
       last_active = activity[-1]
       dist_to_end  = len(scores) - last_active - 3  # the -3 is to account for the size of the activity region (3bp)
       extension = max(0, 2*flank_size - dist_to_end)
   return is_active, extension


def _safediv(a, b):
    """divide while handling zero denominator"""
    try:
        return float(a)/b
    except ZeroDivisionError:
        if a == 0:
            return 0
        raise ValueError('Cannot have a/b with (a == {}, b == {})'.format(a, b))


def _get_active(scores, depth):
   return [_safediv(float(sum(wscore)),sum(wdepth)) > 0.2
           for wscore, wdepth in izip(sliding_window(scores, 3),
                                      sliding_window(depth, 3))]
