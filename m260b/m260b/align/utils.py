"""\
General purpose alignment utilities

"""
from itertools import islice

def _drop(itr, n):
    """\
    drop n elements from the iterator

    """
    for i in xrange(n):
        itr.next()
    return itr

def kmerize(seqitr, size, stride=1, keep_last=False):
    """\
    Turn the base sequence into a sequence of kmers of size `size` with gap `stride`, e.g.

    kmerize('abcdefghijk', 3, 2) -> ['abc', 'cde', 'efg', 'ghi', 'ijk']

    """
    if size <= 0 or stride <= 0:
        raise ValueError('Both size and stride must be > 0')
    seqitr = iter(seqitr)  # allows strings to be passed in directly
    cache_start = -(size-stride) if stride < size else size
    jump_size = stride-size if stride > size else 0
    chunk_size = min(size + stride - size, size) 
    seq = ''.join(islice(seqitr, size))
    while len(seq) == size:
        yield seq
        cache = seq[cache_start:]
        seqitr = _drop(seqitr, jump_size)
        seq = cache + ''.join(islice(seqitr, chunk_size))
        #print('cache: {}  seq:  {}'.format(cache, seq)) 
    if seq and keep_last:
        yield seq  # last one is not full


def compsam(rec1, rec2):
    if rec1.is_unmapped and rec2.is_unmapped:
        return 0
    elif rec1.is_unmapped:
        return 1
    elif rec2.is_unmapped:
        return -1
    if rec1.tid == rec2.tid:
        return -1 if rec1.pos < rec2.pos else 1 if rec2.pos < rec1.pos else 0
    return -1 if rec1.tid < rec2.tid else 1 if rec2.tid < rec1.tid else 0
 
