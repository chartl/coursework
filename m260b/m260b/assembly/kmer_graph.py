"""\
Kmer-graph related code. Sources:

http://www.homolog.us/Tutorials/index.php?p=1.1&s=1
http://www.broadinstitute.org/gatk/guide/tagged?tag=assembly
http://github.com/broadgsa/gatk/   (Haplotype evaluation/HMM; see PairHMM)

"""
from __future__ import print_function

from collections import defaultdict, Counter, namedtuple
from itertools import chain
from operator import itemgetter

import numpy as np

from m260b.align.utils import kmerize
from m260b.assembly.rlm import rlm_readlikelihood
from m260b.traversal.active_window import sliding_window

AssemblyResult = namedtuple('AssemblyResult', ('fail_reason', 'seq', 'event_scores'))
LikCache = namedtuple('LikelihoodCache', ('match', 'insert', 'delete'))
LikPath = namedtuple('LikPath', ('next', 'likelihoods'))


LARGE_VALUE = 2. ** 1020
LARGE_VALUE_LOG = np.log(LARGE_VALUE)

def build_kmergraph(sequences, k=11, min_kmer_count=2, ensure_ref=True):
    """\
    Given a collection of sequences (strings) build a kmer graph of size k.

    For now this can just be a dict.

    """
    graph = defaultdict(Counter)
    counts = Counter()
    sequences = iter(sequences)
    ref = next(sequences)  # special case the first
    graph['$'][ref[:k]] += 1  # sentinel start
    add_sequence(graph, counts, ref, k)
    for seq in iter(sequences):
        add_sequence(graph, counts, seq, k)
    for kmer, count in counts.iteritems():
        # remove all occurrences of it if the count is low
        if count < min_kmer_count and kmer not in {'@', '$'}:
            for weightcount in graph.itervalues():
                weightcount.pop(kmer, None)
            graph.pop(kmer, None)
    if ensure_ref:
        graph['$'][ref[:k]] += 1  # sentinel start
        add_sequence(graph, counts, ref, k)
        graph[ref[-k:]]['@'] += 1  # sentinel end
    return graph


def add_sequence(kgraph, kcounts, seq, k):
    """Add the sequence to the graph"""
    # utility functions for the win
    kto = None
    for kfrom, kto in sliding_window(kmerize(seq, k), 2):
        kgraph[kfrom][kto] += 1
        kcounts[kfrom] += 1
    if kto:
        kcounts[kto] += 1


def build_haplotype(ref, reads, k=11, min_kmer_count=2):
    """\
    Given a reference window and a collection of reads, assemble the
    most likely haplotype, given the reads.

    Note that of all code in this repository, this will be the most
    difficult to generalize to non-haploid organisms.

    """
    graph = build_kmergraph(chain([ref], (read.seq for read in reads)), k, min_kmer_count)
    possible_haplos = _valid_paths(graph)
    haplotype_likelihoods = _evaluate_haplotypes(reads, possible_haplos, del_probs={'open': 0.0005, 'continue': 0.003},
                                                     ins_probs={'open': 0.0005, 'continue': 0.002})
    if len(haplotype_likelihoods) == 0:
        fail_reason = 'no-spanning-reads'
        return AssemblyResult(fail_reason, ref, [])
    best_haplotype = haplotype_likelihoods.keys()[np.argmax(haplotype_likelihoods.values())]
    return AssemblyResult(None, best_haplotype, score_events(haplotype_likelihoods, k))


def _valid_paths(kgraph, node='$', used_nodes=None):
    if node == '$':
        nxt = kgraph[node].keys()[0]  # only one valid first path
        for path in _valid_paths(kgraph, nxt, {nxt}):
            yield nxt[:-1] + path
    elif len(kgraph[node]) == 0:
        yield node[-1]
    else:
        for next_node in kgraph[node]:
            if next_node in used_nodes:
                continue  # no cycles
            for path in _valid_paths(kgraph, next_node, used_nodes | {next_node}):
                if path[-1] == '@':
                    yield node[-1] + path


def _evaluate_haplotypes(reads, putative_haplotypes, del_probs, ins_probs):
    return {haplotype[:-1]: sum(_haplotype_likelihoods(reads, haplotype[:-1],
                                                       del_probs, ins_probs))
            for haplotype in putative_haplotypes}


def _haplotype_likelihoods(reads, putative_haplotype, del_probs, ins_probs):
    for read in reads:
        lik = _c_read_likelihood(read, putative_haplotype) 
        #print(putative_haplotype, read.seq, lik)
        yield lik


_CACHE = None

def _c_read_likelihood(read, haplotype):
    global _CACHE
    if _CACHE is None or _CACHE.match.shape[0] != len(read.seq):
        _CACHE = LikCache(np.zeros((len(read.seq),), dtype=np.float64),
                          np.zeros((len(read.seq),), dtype=np.float64),
                          np.zeros((len(read.seq),), dtype=np.float64))
    read_qual = read.qual or '2' * len(read.seq)
    return rlm_readlikelihood(_CACHE.match, _CACHE.insert, _CACHE.delete,
                              read.seq, haplotype, read_qual,
                              len(haplotype), len(read.seq))


def _read_likelihood(read, haplotype, del_probs, ins_probs):
    """\
    Given a putative haplotype and a read, calculate the likelihood of
    the read coming from the given haplotype. This is very similar to
    smith-waterman except that the calculations are log-probabilities
    as opposed to integer scores; and the goal is to maximize the
    probability.

    `del_probs` and `ins_probs` are dictionaries of
         'open' -> probability of opening such an event
         'continue' -> probability of continuing the event

     ideally we would want
      - maps of context -> floats; these are conditional probabilities
                of an indel event occurring for a given context
                (these would be platform-specific sequence transitions,
                 i.e. TATATA --> TATATATA   0.003; etc)

    :param read: the read
    :param haplotype: the haplotype
    :param del_probs: the probability that the machine skipped a base
    :param ins_probs: the probability that the machine put in some random base

    """
    global _CACHE
    # somehow we need to know whether the last state was an insertion, deletion, or match;
    # and since we're summing we can't just keep track of a 'move' matrix like SW
    if _CACHE is None or _CACHE.match.shape[0] != len(read.seq):
        _CACHE = LikCache(np.zeros((len(read.seq),), dtype=np.float64), 
                          np.zeros((len(read.seq),), dtype=np.float64), 
                          np.zeros((len(read.seq),), dtype=np.float64))
    _init_del = LARGE_VALUE / len(haplotype)
    i2m = 1 - ins_probs['continue']  # mixed terminology here; but these are priors on the transition probabilities
    d2m = 1 - del_probs['continue']
    m2m = 1 - ins_probs['open'] - del_probs['open']
    ins_sum, match_sum = 0., 0.
    for row in xrange(len(haplotype)):
        for col in xrange(len(read.seq)):
            if row is 0:
                lrm, lrd, lri = 0, 0, 0  # can't have anything to do with positions before the start of haplotype
            else:
                lrm, lrd, lri = _CACHE.match[col], _CACHE.delete[col], _CACHE.insert[col]  # previous row, this column
            if col is 0:
                # can't have an insertion or match; only del continuations
                _CACHE.match[col] = _init_del * base_prob(haplotype[row], read, col, transition=d2m)
                _CACHE.delete[col] = lrd * delete_cont_prob(haplotype, row, read, col, del_probs)
                _CACHE.insert[col] = 0
            else:
                _CACHE.match[col] = lrcm * base_prob(haplotype[row], read, col, transition=m2m) + \
                                    lrci * base_prob(haplotype[row], read, col, transition=i2m) + \
                                    lrcd * base_prob(haplotype[row], read, col, transition=d2m)
                _CACHE.delete[col] = lrm * delete_open_prob(haplotype, row, read, col, del_probs) + \
                                     lrd * delete_cont_prob(haplotype, row, read, col, del_probs)
                _CACHE.insert[col] = _CACHE.match[col-1] * insert_open_prob(haplotype, row, read, col, ins_probs) + \
                                     _CACHE.insert[col-1] * insert_cont_prob(haplotype, row, read, col, ins_probs)
            lrcm, lrcd, lrci = lrm, lrd, lri  # previous row, previous column
        ins_sum += _CACHE.insert[-1]
        match_sum += _CACHE.match[-1]
    lik = np.log(ins_sum + match_sum) - LARGE_VALUE_LOG
    #print('{} -> {}'.format(haplotype, lik))
    return lik


def base_prob(hap_base, read, read_offset, transition):
    if hap_base == read.seq[read_offset]:
        # match; probability of outcome is 1-qual
        base_prb = 1. - 10. ** -((ord(read.qual[read_offset])-33)/10.)
    else:
        # mismatch; probability of outcome is qual
        base_prb = 10 ** -((ord(read.qual[read_offset])-33)/10.)
    return transition * base_prb 


def insert_open_prob(hap_seq, seq_offset, read, read_offset, ins_probs):
    """Probability of an insert opening (after a match)"""
    # pretty trivial implementation. Ideally this would be a multi-base thing
    # that knows duplicated sequence is more likely than novel sequence
    # insertion.
    return ins_probs['open'] * 0.25  # 0.25 is base identity


def insert_cont_prob(hap_seq, seq_offset, read, read_offset, ins_probs):
    """Probability of an insertion continuing"""
    return ins_probs['continue'] * 0.25


def delete_open_prob(hap_seq, seq_offset, read, read_offset, del_probs):
    """Probability of a deletion opening (after a match)"""
    return del_probs['open']  # nothing to observe in the read


def delete_cont_prob(hap_seq, seq_offset, read, read_offset, del_probs):
    """Probability of staying in an open deletion state"""
    return del_probs['continue']  # nothing to observe in the read


def _greedy_seq(node):
    return node.most_common(1)[0][0]


def greedy_path(kgraph):
    seq, node, path, pathset, nxt_seq = '$', kgraph['$'], list(), set(), None
    while len(node) > 0:
        path.append(node)
        seq = nxt_seq
        nxt_seq = _greedy_seq(node)
        node = kgraph[nxt_seq]
        if nxt_seq in pathset: 
            # cycle
            return 'cycle', path 
        pathset.add(nxt_seq)
    if '@' in path[-1]:  # made it back to the reference path
        return None, path[:-1]
    return 'no-spanning-reads', path  # didn't assemble through


def path2seq(nodepath):
    seq = _greedy_seq(nodepath[0]) 
    for node in nodepath[1:]:
        seq += _greedy_seq(node)[-1]
    return seq
    path.append(node)
    return path


def score_events(hapliks, k):
    """\
    Now that we have haplotype likelihoods it's time to split them into events. This
    is another graph-construction deal, but this time, instead of counting reads
    that take a transition, we sum the likelihoods (in real space).

    Finally, we take the most likely haplotype and thread it through the graph;
    the score for each element of the haploptype is the phred-likelihood-ratio between
    the branch that includes the most likely haplotype, and the branch that does not.

    Ideally we would marginalize over all haplotypes, but it's not immediately clear
    to me how to demonstrate that two sub-branches line up; e.g.

                  o
                 / \
                /   \
               o     o
              / \    |\
             /   \   o o  <--- same as
 here  -->  o     o  |/
             \   /   o
              \ /   /
               o   /
                \ /
                 o
    """
    _zero = lambda: 0.
    _mkdct = lambda: defaultdict(_zero)
    haplotype_graph = defaultdict(_mkdct)
    best, bestvalue = max(hapliks.items(), key=itemgetter(1))
    for haplotype, likelihood in hapliks.iteritems():
        for kfrom, kto in sliding_window(kmerize(haplotype, k), 2):
            # divide by the largest (which will still be O(10^-20)) to bring things
            # closer to one. The likelihood ratios are unaffected by this.
            haplotype_graph[kfrom][kto] += np.exp(likelihood - bestvalue)
    # now calculate the event likelihoods
    event_liks = -np.ones((1+len(best),), dtype=np.int32)
    offset = k
    for i, (kfrom, kto) in enumerate(sliding_window(kmerize(best, k), 2)):
        outpaths = haplotype_graph[kfrom]
        if len(outpaths) > 1:
            other = sum((value for key, value in outpaths.iteritems() if key != kto))
            lik_float = 10 * (np.log10(outpaths[kto]) - np.log10(other)) + 0.5
            event_liks[k + i] = int(lik_float) if not np.isinf(lik_float) else 5000 * np.sign(lik_float)
            #print('(k + i = {}) {}->({} vs rest)  :: {}'.format(k + i, kfrom, kto, event_liks[k + i]))
    return event_liks
