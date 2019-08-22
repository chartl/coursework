"""\
Tests for kmer-based assembly

"""
import pkg_resources
from pprint import pformat, pprint

import numpy as np
import pysam

from m260b.assembly.kmer_graph import *
from m260b.assembly.kmer_graph import _valid_paths, _evaluate_haplotypes, _read_likelihood, \
                                      _c_read_likelihood
from m260b.io.utils import read_basic_fasta
from m260b.tests.test_pileup import mkread


def test_no_varation_assembly():
   """Test that simple noise is correctly aligned and called"""
   #         0         1         2         3         4         5         6
   #         0123456789012345678901234567890123456789012345678901234567890
   ref =    'ACTTAGAGATCCATATACGATATACTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
   h1 =           'AGATCCATATATGATATACTATATCCCACTCTCGCGTACTACGTTCGATAGCGAT'
   h2 =     'ACTTAGAGATCCATATACGATATATTATATCCCACTCTCGCGTACT'
   h3 =     'ACTTAGAGATCCATATACGATATACTACATCCCACTCTCGCGTACT'
   h4 =                    'TACGATATACTATATCGCACTCTCGCGTACTACGATCGATAGCGAT'
   h5 =                         'TATACTATATCGCACTCTCCCGAACTACGATCGATAGCGAT'
   graph = build_kmergraph([ref, h1, h2, h3, h4, h5], min_kmer_count=0, ensure_ref=False)
   path = greedy_path(graph)[1]
   ref_reconstructed = path2seq(greedy_path(graph)[1]) 
   if ref_reconstructed != ref:
       # convert the graph into a simple dict of dicts
       graph = {k1: {k2: v2 for k2, v2 in v1.iteritems()} for k1, v1 in graph.iteritems()}
       path = [list(a.iteritems()) for a in path]
       raise ValueError('ref != reconstructed:\n{}\n{}'
                        'path:\n{}\ngraph:\n{}'.format(ref, ref_reconstructed, pformat(path), pformat(graph)))

def test_snp_assembly():
    """Test that a simple SNP gets the correct path"""
    ref =  'ACTTAGAGATCCATATACGATATACTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    h1  =  'ACTTAGAGATCCTTATACGATATAGTATATC'
    h2  =   'CTTAGAGATCCTTATACGATATACTATATCC'
    h3  =      'AGAGATCCTTATACGATATACTATATCCCACTCT'
    h4  =           'TCCTTATACGATATACTTTATCCCACTCTCGCGTACT'
    h5  =                             'TATCCCACTCTCGCGTACTAGATCGATAG'
    graph = build_kmergraph([ref, h1, h2, h3, h4, h5], min_kmer_count=0, ensure_ref=False)
    sequence = path2seq(greedy_path(graph)[1])
    expected = 'ACTTAGAGATCCTTATACGATATACTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    if sequence != expected:
        raise ValueError('Did not recover expected path, expected/observed:\n'
                         '{}\n{}'.format(expected, sequence))

def test_del_assembly():
    """Test that a simple deletion gets the correct path"""
    ref =  'ACTTAGAGATCCATATACGATATACTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    h1  =  'ACTTAGAGATCCAT'  +  'AT'
    h2  =   'CTTAGAGATCCAT'  +  'ATA'
    h3  =     'TAGAGATCCAT'  +  'ATACT'
    h4  =         'GATCCAT'  +  'ATACTATATC'
    h5  =             'CAT'  +  'ATACTATATCCCACTCT'
    h6  =               'T'  +  'ATACTATATCCCACTCTC'
    h7  =                                   'ACTCTCGCGTACTACGATCGATAGC'
    graph = build_kmergraph([ref, h1, h2, h3, h4, h5, h6, h7], min_kmer_count=0, ensure_ref=False)
    sequence = path2seq(greedy_path(graph)[1])
    expected = 'ACTTAGAGATCCATATACTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    if sequence != expected:
        raise ValueError('Did not recover expected path, expected/observed:\n'
                         '{}\n{}'.format(expected, sequence))


def test_ins_assembly():
    """Test that a simple insertion gets the correct path"""
    ref =  'ACTTAGAGATCCATATACGATATACTATATCCCAC' + 'TCTCGCGTACTACGATCGATAGCGAT'
    h1  =     'TAGAGATCCATATACGATATACTATAT'
    h2  =            'CCATATACGATATTCTATATCCCAC'
    h3  =                'ATACGATATACTATATCCCACGTG'
    h4  =                      'TATACTATATCCCACGTGGTTCTCGCGT'
    h5  =                                     'GTGGTTCTCGCGTACTACGATCGA'
    graph = build_kmergraph([ref, h1, h2, h3, h4, h5], min_kmer_count=0, ensure_ref=False)
    sequence = path2seq(greedy_path(graph)[1])
    expected = 'ACTTAGAGATCCATATACGATATACTATATCCCACGTGGTTCTCGCGTACTACGATCGATAGCGAT'
    if sequence != expected:
        raise ValueError('Did not recover expected path, expected/observed:\n'
                         '{}\n{}'.format(expected, sequence))


def test_pathological_insertion():
    """Test that we don't die given an insertion that is only partially spanned"""
    ref = 'ACTTAGAGATCCATATACGATATAC'   +   'TATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    h1  = 'ACTTAGAGATGCATATACGATATACTC'
    h2  =    'TAGAGATCCATATACGATATACTCGGC'
    h3  =                            'GGCGGCTTATATCCCACTCTCGCGTAC'
    h4  =                               'GGCTTATATCCCACTCTCGCGTACTACGAT'
    h5  =                                  'TTATATCCCACTCTCGCGTACTACGATCGATAGCGAT'
    graph = build_kmergraph([ref, h1, h2, h3, h4, h5], min_kmer_count=0, ensure_ref=False)
    if greedy_path(graph)[0] is None:
        raise ValueError('Path should not be traversible')


def test_cycles_dont_kill_assembly():
    """\
    Make sure that we don't hit an infinite loop due to any cycles.

    This only works for unix platforms.

    """
    import signal
    def timeout(func, kwargs, duration=3):
        def handler(sig, frm):
            raise ValueError('function {} timed out after {} seconds'.format(func, duration))
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(duration)
        res = func(**kwargs)
        signal.alarm(0)
        return res
    #          0         1         2         3         4         5         6
    #          0123456789012345678901234567890123456789012345678901234567890
    ref =     'ACATCAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTATTATGAC'
    graph = build_kmergraph([ref], k=7, min_kmer_count=0, ensure_ref=False)
    result = timeout(greedy_path, {'kgraph': graph})
    if result[0] is None:
        raise ValueError('Expected non-None result, found {}'.format(result))


def test_path_counting():
   """Check that path counting from the assembly is working properly"""
   #          0    1         2         3         4         5         6         7
   #          567890123456789012345678901234567890123456789012345678901234567890123
   ref =     'TATAGATCTAGCGGCCTATTGCATGTACGTATACGGCAGTCACGTCGTCGCTAATAGCGATCCACTACT'
   h1  =     'TATAGATCTAGCGGCCTATTGCATGTACGT'
   h2  =         'GATCTAGCGGCCTATTGCATGTACGTAT'  +  'TC'
   h3  =             'TAGCCGCCTATTGCATGTACGTAT'  +  'TCACGTCGTCG'
   h4  =                'CGGCCTATTGCATGTACGTAT'  +  'TGACGTCGTCGCTAATAGC'
   h5  =                   'CCTTTTGCATGTACGTAT'  +  'TCACGTCGTCGCTAATAGCGATCC'
   h6  =                             'GTACGTAT'  +  'TCACGTCATCGCTAATAGCGATCCACTACT'
   h7  =                                'CGTAT'  +  'TCACGTCGTCGCTAATAGCGATCCACTACT'
   h8  =                                  'TAT'  +  'TCACGTCGTCGCTAATAGCGATCCACTACT'
   ref_pos = 5
   escores =            '222233333334455555655'  +  '432333322222233333333333222222'
   h1p, h2p, h3p, h4p, h5p, h6p = 5, 9, 13, 16, 19, 29
   reads = [mkread(h, None, '{}M'.format(len(h)), p) # cigar doesn't matter for assembly (!)
            for h, p in [(h1, h1p), (h2, h2p), (h3, h3p), 
                         (h4, h4p), (h5, h5p), (h6, h6p)]]
   q = list(reads[5].qual)
   q[15] = chr(3 + 33)
   reads[5].qual = ''.join(q)
   result = build_haplotype(ref, reads, k=9, min_kmer_count=0)
   expected = 'TATAGATCTAGCGGCCTATTGCATGTACGTAT'  +  'TCACGTCGTCGCTAATAGCGATCCACTACT'
   if result.seq != expected:
       raise ValueError('Observed not equal expected:\n'
                        '{}\n{}'.format(result.seq, expected)) 

def test_valid_path_enumeration():
    """Test that all valid paths can be enumerated"""
    #
    #
    ref = 'TATTCAGATCGCGTATGCAGTATGGGCTATAAACACAGATAGCGAGTCTCAG'
    h1  =     'CAGATCGCGTATGCAGTATGGGCTATATACACAGAT'
    h2  =         'TCGCGTTTGCAGTATGGGCTATAAACACAGATAGCGAG'
    h3  =                      'ATGGGCTATATACACAGATAGGGAGTCTCAG'
    graph = build_kmergraph([ref, h1, h2, h3], 9, min_kmer_count=0)
    pprint(dict(graph))
    expected = {'TATTCAGATCGCGTATGCAGTATGGGCTATAAACACAGATAGGGAGTCTCAG@',
                'TATTCAGATCGCGTATGCAGTATGGGCTATAAACACAGATAGCGAGTCTCAG@',
                'TATTCAGATCGCGTATGCAGTATGGGCTATATACACAGATAGGGAGTCTCAG@',
                'TATTCAGATCGCGTATGCAGTATGGGCTATATACACAGATAGCGAGTCTCAG@'}
    observed = set(_valid_paths(graph))
    if observed != expected: 
        raise ValueError('Sets not equal:\n{}\n{}'.format(observed, expected))


def test_small_snp_assembly():
    """Test a small assembly with a single SNP"""
    ref   =   "TCGGCTATGAAC"
    r1    =   "TCGGCTCT"
    r2    =    "CGGCTCTG"
    r3    =     "GGCTCTGA"
    r4    =      "GCTCTGAA"
    r5    =       "CTCTGAAC"
    def _mk(seq, p):
        return mkread(seq, [20]*len(seq), '{}M'.format(len(seq)), p)
    reads = [_mk(b, a) for a, b in enumerate([r1, r2, r3, r4, r5])]
    graph = build_kmergraph([ref] + [r1, r2, r3, r4, r5], 5, min_kmer_count=0)
    pprint(dict(graph))
    paths = list(_valid_paths(graph))
    expected = {ref + '@', "TCGGCTCTGAAC@"}
    if set(paths) != expected:
        raise ValueError("Path sets not equal: \n{}\n{}".format(sorted(paths), sorted(expected)))
    likelihoods = _evaluate_haplotypes(reads, paths, {'open': 0.005, 'continue': 0.6}, {'open': 0.005, 'continue': 0.6})
    best = likelihoods.keys()[np.argmax(likelihoods.values())]
    if best != "TCGGCTCTGAAC":
        raise ValueError('Best is not better than reference!')


def test_no_passenger_snp():
    """\
    Test for a situation revealed by Issue #11: nearby variation can creep
    in to a non-ref haplotype even if there are only a few variants

    """
    ref_file = pkg_resources.resource_filename('m260b.test_data', 'ref_practice_W_1_chr_1.fasta')
    reads_file = pkg_resources.resource_filename('m260b.test_data', 'practice_w_1.std.bad_region2.bam')
    ref_hdr, reference = read_basic_fasta(ref_file)
    reads = [read for read in pysam.Samfile(reads_file) if not read.is_unmapped]
    haplotype = build_haplotype(reference[1700:1801], reads, k=11, min_kmer_count=2)
    expected = 'CAATATAGCTACGAAACACGTGTGACAACTAAGGCAGGCGATATGAGTCACGGCTTAAATCACTAAGATTGTATTTGGTAAAAGTTGCAGGTGGAGGTCCT'
    if expected != haplotype.seq:
        raise ValueError('Expected haplotype does not match observed:'
                         '\n{}\n{}'.format(expected, haplotype.seq))


def test_read_likelihood():
    """Test that the read likelihood is doing the right thing"""
    ref = 'ACTAGAGATCGTTAGCGTCTCGATCGACGATGACGTTAAGGCCATTAGCGAT'
    # first, a read that matches the ref at Q20 should have a well known score
    r1  =         'TCGTTAGCGTCTCGATCGACGATGACGTTAAGG'
    r2  =         'TCGTTAGCGTCTCTATCGTCGATGTCGTTAAGG'  # snphere
    read1 = mkread(r1, [20]*len(r1), '{}M'.format(len(r1)), 1)
    likelihood = _read_likelihood(read1, ref, {'open': 0.005, 'continue': 0.05},
                                  {'open': 0.005, 'continue': 0.05})
    if likelihood < np.log(0.96 ** len(r1) / len(r1)):
        raise ValueError('Read likelihood is less than it ought to be.'
                         '\nExpected: {}, Observed: {}'.format(np.log(0.98 ** len(r1) / len(r1)),
                          likelihood))
    read2 = mkread(r2, [20]*len(r2), '{}M'.format(len(r2)), 1)
    likelihood2 = _read_likelihood(read2, ref, {'open': 0.005, 'continue': 0.05},
                                  {'open': 0.005, 'continue': 0.05})
    if likelihood2 - likelihood > -0.1:
        raise ValueError('L2={} should be < L1={}'.format(likelihood2, likelihood))


def test_cread_likelihood():
    """test the c read likelihood"""
    def _is_rel_close(a, b, frac):
        return -frac < (a-b)/a < frac
    ref = 'TATTAGA'
    r1  =    'TAG'
    r2  =  'ATCAG'
    r3  = 'TTTAGG'
    for q in xrange(15, 25):
        read = mkread(r1, [q]*len(r1), '{}M'.format(len(r1)), 1)
        r1_lik =_read_likelihood(read, ref, {'open': 0.005, 'continue': 0.05},
                                 {'open': 0.005, 'continue': 0.05})
        r1_clik = _c_read_likelihood(read, ref)
        assert _is_rel_close(r1_lik, r1_clik, 0.01) 
    for q in xrange(15, 25):
        read = mkread(r2, [q]*len(r2), '{}M'.format(len(r2)), 1)
        r2_lik =_read_likelihood(read, ref, {'open': 0.005, 'continue': 0.05},
                                 {'open': 0.005, 'continue': 0.05})
        r2_clik = _c_read_likelihood(read, ref)
        assert _is_rel_close(r2_lik, r2_clik, 0.01)
    for q in xrange(15, 25):
        read = mkread(r3, [q]*len(r3), '{}M'.format(len(r3)), 1)
        r3_lik =_read_likelihood(read, ref, {'open': 0.005, 'continue': 0.05},
                                 {'open': 0.005, 'continue': 0.05})
        r3_clik = _c_read_likelihood(read, ref)
        assert _is_rel_close(r3_lik, r3_clik, 0.05), (r3_lik, r3_clik)
