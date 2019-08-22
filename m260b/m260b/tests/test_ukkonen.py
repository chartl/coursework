from collections import namedtuple
from numpy import array_equal
import numpy as np

from m260b.align.ukkonen import full_sw, banded_sw, _full_sw_matrix, _banded_sw_matrix


def test_full_and_banded_sw():
    """\
    Test that full smith-waterman produces alignments that we would want.

    If anything this has helped to nail down reasonable parameters.

    """
    ref = 'AATTGTACTATTACACTATCGGAGCTAGCTATC'
    x1 =        'ACTATTATACTATC'
    x2 =       'T' +'TATTACACTATCG'
    x3 =    'TTGT' + 'TACACTATCG'
    x4 =         'CTATTATAGGGGGCTATCG'
    x5 =               'CACTATAGGTGCTAGCTATC'  # two snps
    for swname, swfunc in [('full', full_sw), ('banded', banded_sw)]:
	offset, cigar, _, _ = swfunc(ref, x1)
	if offset != 6:
	    raise ValueError('({}): Expected offset 6, found: {}'.format(swname, offset))
	if cigar != '14M':
	    raise ValueError('({}): Expected 14M, found: {}'.format(swname, cigar))
	offset, cigar, _, _ = swfunc(ref, x2)
	if offset != 5:
	    raise ValueError('({}): Expected offset 6, found: {}, {}'.format(swname,offset, cigar))
	if cigar != '1M2D13M':
	    raise ValueError('({}): Expected 1M2D13M, found {}, {}'.format(swname, offset, cigar))
	offset, cigar, _, _ = swfunc(ref, x3)
	if offset != 2:
	    raise ValueError('({}): Expected offset 2, found: {}, {}'.format(swname, offset, cigar))
	if cigar != '4M5D10M':
	    raise ValueError('({}): Expected 4M5D10M, found: {}'.format(swname, cigar))
	offset, cigar, _, _ = swfunc(ref, x4)
	if offset != 7 or cigar != '8M5I6M':
	    raise ValueError('({}): Expected offset 7, 8M5I6M; found: {}, {}'.format(swname, offset, cigar))
	offset, cigar, _, _ = swfunc(ref, x5)
	if offset != 13 or cigar != '{}M'.format(len(x5)):
            full, fm = _full_sw_matrix(ref, x5)
            band, bm = _banded_sw_matrix(ref, x5)
            print(full)
            print(band)
            print(fm)
            print(bm)
	    raise ValueError('({}): Expected offset 13, {}M; found: {}, {}'.format(swname, len(x5), offset, cigar))


def _test_banded_sw_matrix():
   """\
   The banded matrix should produce the exact same matrix as the full sw for large differences

   This test is deprecated and disabled since switching to a fixed c implementation

   """
   ref1 =  'AATGATCTAGAC'
   read1 =  'ATGATCTAG'
   import m260b.align.ukkonen
   m260b.align.ukkonen.BAND_DIFF = 800
   fmat, fmov = _full_sw_matrix(ref1, read1)
   bmat, bmov = _banded_sw_matrix(ref1, read1)
   if not array_equal(fmat, bmat):
       raise ValueError('Banded and full should be the same for a small example.'
                        '\nFull:\n{}\n{}\n\nBanded:\n{}\n{}'.format(fmat, fmov, bmat, bmov)) 


def test_banded_sw():
    """\
    Test that the banded smith waterman captures the same events as the full. These are real examples of screw-ups.
 
    """
    example = namedtuple('TestExample', ('ref', 'read', 'expected_cigar'))
    examples = list()
    ref1 = 'CCCTTTAGTTATGCTTTCTCTTCGGCGGGCGTGGGAC' + 'CGTAATGAGAACTGTACATCAGTCTG'
    read1 =         'TATGCTTTCTCTTCGGCGGGCGTGGGACAAAATCGTAATGAGAACTGTAC'

    cig1 = '28M5I17M'
    examples.append(example(ref1, read1, cig1))
    for example in examples:
        _, cigar, _, _ = banded_sw(example.ref, example.read)
        if cigar != example.expected_cigar:
            import numpy
            numpy.set_printoptions(threshold='nan')
            full_mat, _ = _full_sw_matrix(example.ref, example.read)
            banded_mat, _ = _banded_sw_matrix(example.ref, example.read) 
            raise ValueError('Expected cigar: {}, observed: {}'.format(example.expected_cigar, cigar) + 
                             '\nFull:\n{}\n\nBanded:\n{}'.format(full_mat, banded_mat))


def test_banded_sw_leak():
    """Test that the banded sw implementation is not leaking memory"""
    import os
    import time
    import gc
    from guppy import hpy
    h = hpy()
    proc = os.getpid()
    _Example = namedtuple('TestExample', ('ref', 'read', 'expected_cigar'))
    ref1 = 'CCCTTTAGTTATGCTTTCTCTTCGGCGGGCGTGGGAC' + 'CGTAATGAGAACTGTACATCAGTCTG'
    read1 =         'TATGCTTTCTCTTCGGCGGGCGTGGGACAAAATCGTAATGAGAACTGTAC'

    cig1 = '28M5I17M'
    example = _Example(ref1, read1, cig1)
    def _get_mem():
       size = [t for t in [u.split('\t') for u in open('/proc/{}/status'.format(os.getpid())).read().split('\n')] 
               if t[0] == 'VmSize:']
       print(size)
       return int(size[0][1].split()[0])
    for _ in xrange(5000):
        _ = banded_sw(example.ref, example.read)
    _ = None
    gc.collect()
    heap_init = h.heap()
    initmem = _get_mem()
    for _ in xrange(5000):
        _ = np.zeros((len(example.ref), len(example.read)), dtype=int)
    _ = None
    gc.collect()
    slmem = _get_mem()
    for _ in xrange(5000):
        _ = banded_sw(example.ref, example.read)
    _ = None
    gc.collect()
    curmem = _get_mem()
    if curmem - initmem > 40 and slmem - initmem < 5:
        print('Initial heap:\n{}'.format(heap_init))
        print('-'*40)
        print(h.heap())
        raise ValueError('There is a memory leak. Mem difference: {}kb'.format(curmem - initmem))


def test_shady_haplotype_alignment():
    """Test a real example of a bad alignment"""
    #                                             v                                                          v
    ref = 'AACAACAACAA' +  'CCTGGTCAGGAGTTGAGCCTCCATACTATACTTACTAGTGGTGTACTAACATCCAAACTATTCCCGCGGGACTTAATATGTGATGTCCGCCGTGGTGCGCAATTACGTACGTAGGAAGAGATTGTTATCCAATCTTTTCACGT'
    h1 =  'AACAACAACAACGACAACCTGGTCAGGAGTTGAGCCTCCTTACTATACTTACTAGTGGTGTACTAACATCCAAACTATTCCCGCGGGACTTAATATGTAATGTCCGCCGTGGTGCGCAATTACGTACGTAGGAAGAGATTGTTATCCAATCTTTTCACGT' 
    ref = 'AACAACAAC'  + 'AACCTGGTCAGGAGTTGAGCCTCCATACTATACTT'
    h1 =  'AACAACAACAACGACAACCTGGTCAGGAGTTGAGCCTCCTTACTATACTT' 
    expected_cigar, expected_mismatch = '12M6I32M', 1
    offset, cigar, score, mismatch = banded_sw(ref, h1)
    gain, moves = _banded_sw_matrix(ref, h1)
    assert cigar == expected_cigar, 'E={} != O={}'.format(expected_cigar, cigar)
    assert mismatch == expected_mismatch, 'E={} != O={}'.format(expected_mismatch, mismatch)


def test_shady_haplotype2():
    """Test that a large deletion can actually be deleted"""
    ref = 'CAATCCCCTAGCGGCTCAATCACTGAACCTCCTCCTCTCCGGGGCGTTGGCGTCTTCTTTTATGTGAGAAGAATAATTACCCCTAGCGGCGTTAACAGTTGGGTG'
    h1  = 'CAATCCCCTAGCGGC'                                 +                                       'GTTAACAGTTGGGTG'
    expected_cigar = '15M75D15M'
    offset, cigar, score, mismatch = banded_sw(ref, h1, not_in_ref_penalty=40)
    foffset, fcigar, fscore, fmismatch = full_sw(ref, h1, lenient=True)
    assert fcigar == expected_cigar
