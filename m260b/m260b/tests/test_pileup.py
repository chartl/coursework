"""\

Make sure the pileup objects are doing the right thing

"""
import re

from pysam import AlignedRead

import m260b.traversal.pileup
from m260b.traversal.pileup import PileupElement, Pileup, make_pileups

_READS_MADE = 0
def mkread(seq, quals, cigar, pos=None):
    global _READS_MADE
    r = AlignedRead()
    r.qname = 'testread{}'.format(_READS_MADE)
    r.seq = seq
    quals = [22] * len(seq) if quals is None else quals
    r.qual = ''.join([chr(q+33) for q in quals]) 
    r.cigarstring = cigar
    if pos:
        r.pos = pos
    _READS_MADE += 1
    return r

def test_sequence_access_is_correct():
    base_seq = 'ACTTAGTATCCTATTGACTATGGATCTAC'
    for test_cigar in ['29M', '20M4D9M', '15MD59M', '1M16D28M', '3M4D8M5D12M2D9M']:
        element = PileupElement(mkread(base_seq, None, test_cigar), 0)
        bases = []
        while not element.is_last_position():
            bases.append(element.base())
            element._next_position()
        bases.append(element.base())
        bases = filter(lambda m: m != None, bases)
        if ''.join(bases) != base_seq:
            raise ValueError('test_sequence_access_is_correct failed: {} != {}'.format(''.join(bases), base_seq))


def test_cigar_update_properly_updates_cigar():
    base_seq = 'ACTTAGTATCCTATTGACTATGGATCTAC'
    for offset in [2, 3, 6, 8, 12, 13, 18]:
        cigar = '{}M{}D{}M'.format(offset-1,6,29-offset+1)
        element = PileupElement(mkread(base_seq, None, cigar), 0)
        for _ in xrange(offset):
            element._next_position()
        if element._cigarelement()[0] != 2:
            raise ValueError('Not in deletion element at offset {} of read with cigar {}'.format(offset, cigar))
        cigar = '{}M{}I{}M'.format(offset-1, 3, 26-offset+1)
        element = PileupElement(mkread(base_seq, None, cigar), 0)
        for _ in xrange(offset):
            element._next_position()
        if element._cigarelement()[0] == 1:
            raise ValueError('Did not properly skip insertion element at offset {} of read with cigar {}'.format(offset, cigar))


def test_deletions_dont_return_bases():
    base_seq = 'ACTTAGTATCCTATTGACTATGGATCTAC'
    for offset in [2, 3, 5, 8, 13, 19]:
        for length in [1, 2, 5, 8]:
            cigar = '{}M{}D{}M'.format(offset - 1, length, 29 - offset + 1)
            element = PileupElement(mkread(base_seq, None, cigar), 0)
            bases = []
            while not element.is_last_position():
                bases.append(element.base())
                element._next_position()
            bases.append(element.base())
            none_elem = filter(lambda t: t == None, bases)
            if len(none_elem) != length:
                raise ValueError('pileup elements not returning None for within-deletion bases:\n{}\n{}'.format(cigar, bases))


def test_start_at_offset():
    """\
    Test that starting a cigar element at an offset doesn't corrupt the base sequence. Deletions are tricky.

    Ref:   ACTGTTATCCCATAATCGAGAGAGCTATTAGCGATCCCTATCTAGAGGCCGGAAGGCCATATTAGAC
    Read:              TAATCGAGAGAGCTATT-----TCCCTAT
    Position:                             ^

    You have to know that that offset (position - read_start) into the read is in the middle of that deletion

    Similarly:
 
    Ref:   ACTGTTATCCCATAAT----CGAGAGAGCTATTAGCGATCCCTATCTAGAGGCCGGAAGGCCATATTAGAC
    Read:         TCCCATAATCGAGCGAGAGAGCTATTAGCGATCCC
    Position:                  ^

    You have to know that offset (position - read_start) really means the reference shift, and so takes you
    *beyond* the insertion sequence into the read

    """

    base_seq = 'ACTTAGTATCCTATTGACTATGGATCTAC'
    for ref_delta in xrange(23):
        for cigar in ['29M', '5M4D24M', '14M5D15M', '10M4I15M', '10M2D3M2I14M']:
            element = PileupElement(mkread(base_seq, None, cigar), ref_delta)
            cigars = [[int(t[0]), t[1]] for t in re.findall('([0-9]+)([MDISHX])', cigar)]
            expected_seq = []
            base_counter = 0
            seq_offset = 0
            # ref delta means at the Kth *reference* position from the start of the read
            #    insertion sequences occuring before K do not count, while sequences following
            #    should also be skipped.
            #    deletion sequences before K *do* count, and sequences following will return None
            while base_counter < ref_delta:
                cigars[0][0] -= 1
                if cigars[0][1] != 'I':
                    base_counter += 1
                if cigars[0][1] != 'D':
                    seq_offset += 1
                if cigars[0][0] == 0:
                    cigars.pop(0)
            while seq_offset < len(base_seq):
                cigars[0][0] -= 1
                if cigars[0][1] == 'I':
                    seq_offset += 1
                elif cigars[0][1] == 'D':
                    expected_seq.append(None)
                else:
                    expected_seq.append(base_seq[seq_offset])
                    seq_offset += 1
                if cigars[0][0] == 0:
                    cigars.pop(0)
            seq = []
            while not element.is_last_position():
                seq.append(element.base())
                element._next_position()
            seq.append(element.base())
            def n2d(s):
                return s or 'D'
            if ''.join(map(n2d, expected_seq)) != ''.join(map(n2d, seq)): 
                raise ValueError("Sequences do not match.\n"
                                 "Expected: {}\nObserved: {}"
                                 "\nRead_seq:{}\nCigar:{}\nDelta:{}".format(expected_seq, seq, 
                                                                               base_seq, cigar,ref_delta))


def test_base_counts():
    """\
    Test that the pileup object and iterator are doing things by testing the base count functionality.

    We create the following situation:

              0         1         2         3         4
    Pos:      012345678901234567890123456789012345678901
    Ref:      ACTGGGATCTCTTCGACTATATCCGAGCGATCGATCGATGGT
    Read1:                   ACTATATCCGAGCGAT
    Read2:                     TATATCCGAGGGATCG
    Read3:                     TATAT---AGGGATCGATC
    Read4:                      ATAT---AGCGATCGACCG
    Read5:                       TATCCGAGGXATCG

    Where the 'X' in read 5 means GGG, where 2 Gs have been inserted with respect to the reference.

    This leads to counts of
    15 - A:1
    16 - C:1
    17 - T:3
    18 - A:4
    19 - T:5
    20 - A:5
    21 - T:5
    22 - C:3, None:2
    23 - C:3, None:2
    24 - G:3, None:2
    25 - A:5
    26 - G:5
    27 - C:2, G:3
    28 - G:5
    29 - A:5
    30 - T:5
    31 - C:4
    32 - G:4
    33 - A:2
    34 - C:1, T:1
    35 - C:2
    36 - G:1

    """
    reference = 'ACTGGGATCTCTTCGACTATATCCGAGCGATCGATCGATGGT'
    rpos = [0] * 5
    r1seq, r1cig, rpos[0] = 'ACTATATCCGAGCGAT', '16M', 15 
    r2seq, r2cig, rpos[1] = 'TATATCCGAGGGATCG', '16M', 17
    r3seq, r3cig, rpos[2] = 'TATATAGGGATCGATC', '5M3D11M', 17
    r4seq, r4cig, rpos[3] = 'ATATAGCGATCGACCG', '4M3D12M', 18
    r5seq, r5cig, rpos[4] = 'TATCCGAGGGGGATCG', '10M2I4M', 19
    reads = [mkread(r1seq, None, r1cig), mkread(r2seq, None, r2cig), mkread(r3seq, None, r3cig),
             mkread(r4seq, None, r4cig), mkread(r5seq, None, r5cig)]
    for pos, read in zip(rpos, reads):
        read.pos = pos
    expected_counts = {15: {'A': 1}, 16: {'C': 1}, 17: {'T': 3}, 18: {'A': 4}, 19: {'T': 5},
                       20: {'A': 5}, 21: {'T': 5}, 22: {'C': 3, None: 2}, 23: {'C': 3, None: 2},
                       24: {'G': 3, None: 2}, 25: {'A': 5}, 26: {'G': 5}, 27: {'C': 2, 'G': 3},
                       28: {'G': 5}, 29: {'A': 5}, 30: {'T': 5}, 31: {'C': 4}, 32: {'G': 4},
                       33: {'A': 2}, 34: {'C': 1, 'T': 1}, 35: {'C': 2}, 36: {'G': 1}} 
    for pileup in make_pileups(reads, reference, offset=0, contig='blah'):
        expected = expected_counts.get(pileup.position, dict())
        if len(expected) == 0:
            if len(pileup) > 0:
                raise ValueError('Expected 0 counts at {} but pileup size was {}'.format(pileup.position, len(pileup)))
            continue
        for base, count in expected.iteritems():
            if pileup.base_counts[base] != count:
                raise ValueError('Count incorrect for offset '
                                 '{}: {}:{} != {}:{}'.format(pileup.position, base, count,
                                                             base, pileup.base_counts[base]))
