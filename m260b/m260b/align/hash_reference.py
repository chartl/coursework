"""\
Alignment utilities for creating hashes of a fasta reference file

"""
import cPickle

import mmh3  # hashing library

from m260b.io import utils as ioutils
from m260b.align import utils as alnutils


def build_hashtable(contig_seq, kmer_size, stride):
    """\
    Build a hashtable (mapping from an integer to an offset) for the given
    contig sequence. The hashtable consists of a low-overhead unique seed map,
    a higher-overhead collision table, and a metadata table giving information
    about the table

    """
    singlehit = dict()
    collision = dict()
    for idx, kmer in enumerate(alnutils.kmerize(contig_seq, kmer_size, stride)):
        hash = mmh3.hash(kmer)
        offset = stride * idx
        if hash in collision:
            collision[hash].append(offset)
        elif hash in singlehit:
            collision[hash] = [singlehit[hash], offset]
        else:
            singlehit[hash] = offset
    # prune out the collisions
    for hash in collision:
        del singlehit[hash]
    info = {'kmer_size': kmer_size, 'stride': stride, 'unambiguous': len(singlehit),
            'ambiguous': sum((len(v) for v in collision.values()))}
    return singlehit, collision, info


def verify_hash(contig_str, offset, size, expected_value):
    """\
    Verify that a previously computed hash based on a given offset is as
    expected

    """
    return mmh3.hash(contig_str[offset:(offset + size)]) == expected_value


def save_hash(single, collision, info, file):
    pickler = cPickle.Pickler(open(file, 'wb')) 
    pickler.dump([single, collision, info])

def load_hash(hashfile):
    reflist = cPickle.Unpickler(open(hashfile, 'rb')).load()
    return (reflist[0], reflist[1], reflist[2])
