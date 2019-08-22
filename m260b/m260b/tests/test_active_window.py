"""\
Tests for active window functionality

"""
from collections import namedtuple
import pkg_resources
import pysam

import numpy as np

from m260b.io.utils import read_basic_fasta
from m260b.tests.test_pileup import mkread
from m260b.traversal.active_window import get_active_score, _get_active, get_activity, \
                                          _advance_window, _update_scores, _advance_scores, \
                                          active_regions 


def test_active_score():
    """\
    Test that the active score function is working properly. Set up some situations:
             0            
             0
             0  1         2         3         4         5         6         7
             789012345678901234567890123456789012345678901234567890123456789012|
    Ref:     CCTATGTGTGTGCTATACGATAGCGATCCAGTACGGTAAACGTATTCCTGATGGAACTGCGCTCTCGAGATCAGTCATG
    Read1:        GTGTGTGCTATATGATAGCGATCCACTAC
    Exp:          00000000000010000000000001000
    Read2:                       TAGCGAT-----ACGGTAAACGTATTCCTGAT
    Exp:                         0000003000030000000000000000000
    Read3:                                           CGTATTCCXGATGGAACTGCGCTC  (X = 'TCCCCC')
    Exp:                                             0000000030000000000000000
    Read4:      ATGTCTGTGCT------TAGCGATCCAGT
    Exp:        00001000003000003000000000000
    Read5:                                          (X = 'GTG')    AACTXCGCTCACGAGATCAGTCT
    Exp:                                                           000040000010
    Read6: TACCTATGTGTGTGCTTTACGATAGC
    Exp:     000000000000001000000000

    """
    Example = namedtuple('Example', ['seq', 'cigar', 'pos', 'expected'])
    offset = 7
    ref = 'AGTACAA' + 'CCTATGTGTGTGCTATACGATAGCGATCCAGTACGGTAAACGTATTCCTGATGGAACTGCGCTCTCGAGATCAGTCATG'
    examples = [Example('GTGTGTGCTATATGATAGCGATCCACTAC', '29M', 12, map(int, '00000000000010000000000001000')),
           Example('TAGCGATACGGTAAACGTATTCCTGAT', '7M5D20M', 27, map(int, '00000030000300000000000000000000')),
           Example('CGTATTCCTCCCCCGATGGAACTGCGCTC', '9M5I15M', 47, map(int, '000000004000000000000000')),
           Example('ATGTCTGTGCTTAGCGATCCAGT', '11M6D12M', 10, map(int, '00001000003000003000000000000')),
           Example('AACTGTGCGCTCACGAGATCAGTCT', '5M2I18M', 61, map(int, '000040000010')),
           Example('TACCTATGTGTGTGCTTTACGATAGC', '26M', 5, map(int, '000000000000001000000000'))]
    for i, example in enumerate(examples):
        read = mkread(seq=example.seq, cigar=example.cigar, pos=example.pos, quals=None)
        read_offset = read.pos - offset
        active_score = get_active_score(ref, read, offset, 3192).tolist()
        if not all([a == b for a, b in zip(active_score, example.expected)]):
             raise ValueError('Error in example {}. Observed != expected: \n'
                      '{}\n{}'.format(1+i, ''.join(map(str, active_score)),
                                      ''.join(map(str, example.expected))))


def test_get_activity():
    """\
    Test that get_activity performs the sums correctly

    """
    # basically create some depth and score arrays
    Example = namedtuple('Example', ('score', 'depth', 'active_expected', 'shift_expected'))
    examples = [
        Example([0, 0, 0, 0, 1, 3, 1, 0, 0], [5, 5, 8, 8, 12, 14, 15, 15, 15], [False] * 7, 0),
        Example([0, 0, 3, 5, 0, 0, 0, 9, 0], [5, 7, 8, 10, 8, 7, 8, 8, 8], [False, True, True, False, False, True, True], 8),
        Example([31, 30, 30, 0, 0, 0, 0, 0, 0], [120, 180, 150, 150, 130, 80, 40, 20, 40], [True] + [False] * 6, 2),
        Example([0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [False] * 7, 0)]
    for i, example in enumerate(examples):
        activity = _get_active(np.array(example.score), np.array(example.depth))
        active, shift = get_activity(np.array(example.score), np.array(example.depth), 4)
        if not all([a == b for a, b in zip(activity, example.active_expected)]):
            raise ValueError('Error in activity calculation for example {}:\n'
                             'Observed != Expected:\n{}\n{}'.format(1+i, activity, example.active_expected))
        if shift != example.shift_expected:
            raise ValueError('Error in window shift calculation for example {}:\n'
                             '{} != {}'.format(1+i, shift, example.shift_expected))
        if active != any(activity):
            raise ValueError('How did you screw this up, man?')


def test_advance_window():
    """\
    Test that advance window is doing the right thing

    """
    offset = 21273
    wend = 21294
    pool = list()
    for start_pos in xrange(21250, 21274):
        for seqsize in [20, 25, 30, 35, 40, 45, 50]:
            if start_pos + seqsize >= 21273:
                pool.append(mkread(seq='A'*seqsize, cigar='{}M'.format(seqsize), pos=start_pos, quals=None))
            for dsize in [1, 3, 4, 5, 8]:
                if start_pos + seqsize + dsize >= 21273:
                    pool.append(mkread(seq='A'*seqsize, cigar='{}M{}D{}M'.format(10, dsize, seqsize-10), quals=None))
            for isize in [1, 2, 5, 10]:
                if start_pos + seqsize - isize >= 21273:
                    pool.append(mkread(seq='A'*seqsize, cigar='{}M{}I{}M'.format(5, isize, seqsize-5-isize), quals=None))
    pool_copy = [read for read in pool]
    old_offset = offset
    pool, offset, new_wend = _advance_window(pool, offset, wend, 10, 23130)
    # make sure everything left overlaps the new offset
    if not all([read.pos + read.alen >= offset for read in pool]):
        raise ValueError('Error in pool advancement: some reads do not overlap the offset:\n'
                         '{}'.format([read for read in pool if read.pos + read.alen < offset]))
    # make sure everything not in the pool does not overlap the new offset
    kept = {(read.pos, read.alen, read.cigarstring) for read in pool}
    dropped = [read for read in pool_copy if (read.pos, read.alen, read.cigarstring) not in kept]
    if not all([read.pos + read.alen <= offset for read in dropped]):
       raise ValueError('Error in pool advancement; some reads which overlap the offset fail'
                        ' to make it into the pool:\n'
                        '{}'.format('\n'.join([str(read) for read in dropped if read.pos + read.alen >= offset])))
    if offset != old_offset + 11:
        raise ValueError('Error in pool advancement: the new offset does not equal the old offset but shifted by flank+1')

    # now check the edge case
    pool = [read for read in pool_copy]
    pool, offset, new_wend = _advance_window(pool, offset, wend, 10, 21303)
    if new_wend != 21302:
        raise ValueError('The window advance does not properly truncate window size at end of contig')


def test_advance_scores():
    """\
    Test that the score advancement is functioning properly. Set up the following situation:

    scores:          000000112233333444443
    window:        12|---------o---------|32
    reads:                 ------------->
                             ------------>
                               -------------->
                                    ------------->
    next window:             22|---------o---------|42
    next scores:               333334444432222111100

    the reads are all 100% mismatch so that we can validate using sums

    """
    offset, win_end = 12, 32
    ref = 'A'*60
    r1, r1o = 'T'*14, 6
    r2, r2o = 'T'*13, 8
    r3, r3o = 'T'*15, 10
    r4, r4o = 'T'*14, 15
    pool = [mkread(r, None, '{}M'.format(len(r)), offset + ro)
             for r, ro in zip([r1, r2, r3, r4], [r1o, r2o, r3o, r4o])]
    refscore = np.array(map(int, '000000112233333444443'))
    depth = np.array([r for r in refscore])
    expected = map(int, '333334444432222111100')
    refscore, depth = _advance_scores(pool, 22, 42, 10, ref, refscore, depth)
    if not all([a == b for a, b in zip(refscore.tolist(), expected)]):
        raise ValueError('Observed != Expected\n{}\n{}'.format(refscore.tolist(), expected))


def test_advance_scores2():
    """\
    A second test for the advancement of the window

 
    scores:        00001112223339344333322
    window:      23|----------o----------|45
    reads:             ------------>  (27, 0, 13)
                          ------|------>    (30, 3, 14)
                             ---|-------->    (33, 3, 13)
                                  -------------->   (38, 0, 15)  
    next_window:            34|----------o----------|56
    next_scores:              33934433332211111110000

    """
    offset, win_end = 23, 45
    ref = 'A' * 60
    read1 = mkread('T'*13, None, '13M', 27)
    read2 = mkread('T'*17, None, '7M3I7M', 30)
    read3 = mkread('T'*16, None, '4M3I9M', 33)
    read4 = mkread('T'*15, None, '15M', 38)
    pool = [read1, read2, read3, read4]
    refscore = np.array(map(int, '00001112223339344333322'))
    depth = np.array(map(int, '00001112223333344333322'))
    expected = map(int, '33934433332211111110000')
    refscore, depth = _advance_scores(pool=pool, offset=34, win_stop=56, 
                                      flank=11, reference=ref, 
                                      refscore=refscore, depth=depth)
    if not all([a == b for a, b in zip(refscore.tolist(), expected)]):
        raise ValueError('Observed != Expected\n{}\n{}'.format(refscore.tolist(), expected))


def test_advance_scores_reads_bigger_than_window():
    """\
    Advance scores broke when reads were larger than the window. This test
    ensures we shall never regress.

    """
    ref = 'GGTGGAGTAGCCCACGTAGATGCGACCATCC' + 'T' * 100  # suffix doesn't matter
    offset, win_stop, flank = 0, 20, 10
    read1 = mkread('GATGGAGTAGCCCACGTAGATGCGACCATCCTCTCCGGGCTTAGCGGTCT',
                    None, '50M', 0)
    read2 = mkread('CCCACGTAGATGCGACCATCCTCTCCGGGCTTAGCGGTCTTTCTATTCAT',
                   None, '50M', 10)
    read3 = mkread('ACGTAGATGCGACGATCCTCTCCGGGCTTAGCAGTCTTTCTATTCATGTG', # note snp: GATC not CATC
                   None, '50M', 13)
    read4 = mkread('ACGTAGATGCGACCATCCTCTCCGGGCTTAGCGGTCTTTCTATTCATGTG',
                   None, '50M', 13)
    read5 = mkread('CGTAGATGCGACCATCCTCTGCGGGCTTAGCGGTCTTTCTATTCATGTGA',
                   None, '50M', 14)
    read6 = mkread('TAGATGCGACCATCCTCTCCGGGCTTAGCGGTCTTTCTATTCATGTGAGG',
                   None, '50M', 16)
    refscore = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    expected = [0] * 21
    expected[-5] = 1
    refscore, depth = _advance_scores(pool=[read1, read2, read3, read4, read5, read6],
                                      offset=10, win_stop=30, flank=10, reference=ref,
                                      refscore=refscore, depth=np.array([4] * 21))
    if not all([a == b for a, b in zip(refscore.tolist(), expected)]):
        raise ValueError('Observed != Expected\n{}\n{}'.format(refscore.tolist(), expected))


def test_advance_scores_nasty_del_softclip():
    """Test a situation that truly failed: deletion followed by softclip (not really allowed but...)"""
    read1 = mkread('AAAGCGGTTCACAAGACGCCGGACGTATGAGTTGAGAG'+'CTATAAAGTAAA', None, '38M3D10M2S', 2)
    read2 = mkread('AAAGCGGTTCACAAGACGCCGGACGTATGAGTTGAGAG'+'CTATAAAGTAAA', None, '38M3D12M', 2)
    ref =        'TCAAAGCGGTTCACAAGACGCCTGACGTATGAGTTGAGTGGAACGATTTAGTATCATATCTTGGGACGGTCAAATAGACTGTACCCTTCC'
    offset = 0
    win_end  = 50
    score_start = 24
    refscore = np.array(map(int, ('0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 1 0 0 1 2 0 0'
                                  '0 3 0 2 2 0 0 0 0 0 0 0 0 0').split()), dtype=int)
    depth = np.array(map(int, ('1 17 17 17 16 15 14 14 14 14 14 14 13 13 13 12 10  9  8  7  7  7  7  7  7'
                               '0  7  7  7  6  6  6  5  4  4  4  3  3  3  3  2  2  1  0  0  0  0  0  0  0'
                               '0').split()), dtype=int) 
    scores1 = get_active_score(ref, read1, offset, win_end, score_start)
    scores2 = get_active_score(ref, read2, offset, win_end, score_start)
    if scores1.shape[0] != scores2.shape[0] or not all([a == b for a, b in zip(scores1, scores2)]):
        raise ValueError('Different shapes for 10M2S and 12M'
                         '\n{}\n{}'.format(scores1, scores2))
    _advance_scores(pool=[read1], offset=25, win_stop=win_end, flank=25, reference=ref, refscore=refscore, depth=depth) 


def test_update_scores():
    """Test that the code to update a ref score array is working properly"""
    ref = np.zeros((36,), dtype=int)
    depth = np.array([_ for _ in ref])
    ascore = np.array([3, 1, 2, 0, 0, 6])
    for offset in xrange(31):
        pref = [_ for _ in ref]
        prev = ref[offset + 5]
        pdep = depth[offset + 5]
        _update_scores(ref, depth, ascore, offset)
        if ref[offset + 5] != 6 + prev:
            raise ValueError('Error updating scores at offset {}\n'
                             'prev: {}\n'
                             'next: {}'.format(offset, pref, ref.tolist()))
        assert depth[offset + 5] == 1 + pdep


def test_bad_region():
    """Bad region identified in issue #8 should actually trigger"""
    ref_file = pkg_resources.resource_filename('m260b.test_data', 'ref_practice_W_1_chr_1.fasta')
    read_file = pkg_resources.resource_filename('m260b.test_data', 'practice_w_1.std.bad_region1.bam')
    ref_hdr, reference = read_basic_fasta(ref_file) 
    read_iter = pysam.Samfile(read_file)
    chr = ref_hdr[1:].strip()
    areg = list(active_regions(read_iter, reference, chr, start_offset=0, flank=30, dfrac=1.0))
    found = False
    for region, reads in areg:
        found |= region.start <= 5769 <= region.stop
    if not found:
        raise ValueError('Window did not open around variant')
