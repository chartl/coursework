"""\
Utility functions for hashing reads (and aligning them to a reference)

"""
from collections import Counter, namedtuple
import mmh3 
from itertools import izip
import sys

from pysam import AlignedRead

from m260b.align.utils import kmerize
from m260b.align.ukkonen import full_sw, banded_sw
import m260b.debug.debug
from m260b.debug.debug import debug


OffsetInfo = namedtuple('OffsetInfo', ('kmer_offset', 'ref_offset', 'implied_start'))
HashInfo = namedtuple('HashInfo', ('seq', 'fwd', 'rev'))
AlignmentInfo = namedtuple('AlignmentInfo', ('offset', 'cigar', 'reversed', 'mismatches'))
AlignmentScore = namedtuple('AlignmentScore', ('start_offset', 'score', 'cigar', 'mm')) 


USE_SW = True 
MIN_FASTPATH_HITS=20


def _aln_info(kmer_offset, ref_offset):
    return OffsetInfo(kmer_offset, ref_offset, ref_offset-kmer_offset) if ref_offset else None


@debug
def hash_transform(index, paired_reads):
    """\
    Given a reference hash index `index` and a pair of reads, return the read pairs
    along with their hash transform, that is, the list of mapping locations of each of the
    kmers in each read.

    """
    singlehit, collision, info = index
    read1, read2 = paired_reads
    kmer_size = info['kmer_size']
    def transform(read):
        return [_aln_info(idx, singlehit.get(mmh3.hash(kmer), None)) 
                for idx, kmer in enumerate(kmerize(read, kmer_size, stride=1))]
         # note - homework 1 calls for reversing and not reverse complementing. Odd.
    return (HashInfo(read1, transform(str(read1.seq)), 
                     transform(str(read1.seq)[::-1])),
            HashInfo(read2, transform(str(read2.seq)), 
                     transform(str(read2.seq)[::-1])))  # oh god it looks like LISP

@debug
def score_read_mismatch(ref, read, pos):
    """\
    Score the read against the reference. For now just use the # of mismatches.
    In the future switch over to Ukkonen.

    """
    mismatches = (a != b for a, b in izip(ref[pos:(pos + len(read))], read))
    mismatches = [idx for idx, mismatch in enumerate(mismatches) if mismatch]
    return AlignmentScore(pos, len(mismatches), '{}M'.format(len(read)), len(mismatches))


@debug
def score_read(ref, read, pos, use_sw=True, sw_alg=banded_sw):
    if use_sw and USE_SW:
        ref_window = ref[max(pos-20, 0):min(pos+len(read)+20, len(ref)-1)]
        window_offset, cigar, score, mismatches = sw_alg(ref_window, read.seq)
        if window_offset is None:
            return score_read_mismatch(ref, read, pos)
        ref_delta = pos - max(pos-20, 0)
        return AlignmentScore(pos - ref_delta + window_offset, -score, cigar, mismatches)
    return score_read_mismatch(ref, read, pos) 

@debug
def best_alignment_single(ref_seq, index, read_seq, read_hash, max_edit, max_indels, min_seeds, max_hits):
    """Determine the best alignment for a single read and orientation See `best_alignment`."""
    if read_seq == '':
        return None
    start_pos = Counter((info.implied_start for info in read_hash if info))
    if len(start_pos) == 0 or len(start_pos) > max_hits:
        return None
    if max(start_pos.values()) < min_seeds:
        return None
    if len(start_pos) == 1 and start_pos.values()[0] > MIN_FASTPATH_HITS:
        # only one hit, multiple hashes, this is just straightforward. Count mismatches.
        best = score_read(ref_seq, read_seq, start_pos.keys()[0], use_sw=False)
    elif len(start_pos) == 1:
        # one hit but only a few counts, could be noisy. Use SW here.
        best = score_read(ref_seq, read_seq, start_pos.keys()[0], use_sw=True)
    elif len(start_pos) == 2:
        # two hits implies there's an indel, so just run SW on the smaller position
        small_start = min(start_pos.keys())
        best = score_read(ref_seq, read_seq, small_start, use_sw=True)
    else:
        #print('multiple hits found:\n{}'.format(start_pos))
        scores = [score_read(ref_seq, read_seq, pos) for pos in start_pos]  # note: probably only need to check best 2
        best = min(scores, key=lambda r: r.score)
    return best if best.score < max_edit else None
     

@debug
def best_alignment(ref, index, hash1, hash2, max_edit, max_indels, min_seeds, max_hits):
    """\
    Determine the best alignment from the hash-transformed reads. See `align_pairs`.
    
    :param index: the reference index
    :param hash1: the HashInfo for the first read of the pair
    :param hash2: the HashInfo for the second read of the pair
    :param max_edit: maximal edit distance
    :param max_indels: maximal indels
    :param min_seeds: minimum hits
    
    :return: the best alignment (AlignmentInfo) for each read

    Basic cases:
      - All implied offsets the same
      - Implied offsets close together but different
      - Differing implied offsets

    Advanced cases:
      - Any basic + high mismatch rate

    """
    # the structure of this function is because in general we can use the alignment of
    # one read to inform the alignment of its mate. For now, ignore the information
    # that they are paired and just consider them separately.

    # TODO eventually kill off the [::-1] in favor of reverse complement, but HW1 requires only reverse
    r1_fwd, r1_rev, r2_fwd, r2_rev = None, None, None, None
    if hash1.seq.seq:
        r1_fwd = best_alignment_single(ref, index, hash1.seq, hash1.fwd, max_edit, max_indels, min_seeds, max_hits)
        r1_rev = best_alignment_single(ref, index, hash1.seq[::-1], hash1.rev, max_edit, max_indels, min_seeds, max_hits)
    if hash2.seq.seq:
        r2_fwd = best_alignment_single(ref, index, hash2.seq, hash1.fwd, max_edit, max_indels, min_seeds, max_hits)
        r2_rev = best_alignment_single(ref, index, hash2.seq[::-1], hash2.rev, max_edit, max_indels, min_seeds, max_hits)
    def get_aln_info(fwd, rev, size, ref_end):
        if fwd and rev:
            aln = AlignmentInfo(offset = fwd.start_offset if fwd.score > rev.score else rev.start_offset,
                                   reversed=rev.score >= fwd.score, cigar=rev.cigar if rev.score > fwd.score else fwd.cigar, 
                                   mismatches=fwd.mm if fwd.score > rev.score else rev.mm)
        elif fwd:
            aln = AlignmentInfo(offset=fwd.start_offset, reversed=False, cigar=fwd.cigar,
                                   mismatches=fwd.mm)
        elif rev:
            aln = AlignmentInfo(offset= rev.start_offset, reversed=True, cigar=rev.cigar,
                                   mismatches=rev.mm)
        else:
            aln = None
        if aln and (aln.offset + size >= ref_end or aln.offset < 0):
            aln = None
        if aln:
            cigarcount = Counter(aln.cigar)
            if cigarcount['I'] + cigarcount['D'] > max_indels:
                aln = None
        return aln
    r1_aln = get_aln_info(r1_fwd, r1_rev, len(hash1.seq), len(ref))
    r2_aln = get_aln_info(r2_fwd, r2_rev, len(hash1.seq), len(ref))
    return r1_aln, r2_aln


@debug
def align_pairs(ref, index, paired_reads, read_group, max_edit_dist=5, max_indels=2, min_seeds=3, max_hits=8):
    """\
    Given a reference hash index `index` and a sequence of paired reads, return a sequence of
    SAMRecords which contain the information of the best alignment of the paired-end reads
    according to our hash.

    The `best` alignment is the alignment which minimizes the edit distance from the reference.
    However, certain parameters determine whether a read is considered for alignment or not:

    :param max_edit_dist: The maximum edit distance from the reference. If the estimated
      edit distance (from the hashing stage) is > `max_edit_dist`, ignore the read.
    :param max_indels: The maximum number of indels to consider. If the seeds imply that
      there are more than `max_indels` indel events in the read, ignore the read.
    :param min_seeds: The minimal number of seeds to consider. If a read hash contains
      < `min_seeds` hits to the index, ignore the read.

    :param index: the reference index
    :param paired_reads: a sequence of paired reads (Bio.SeqRecord.SeqRecord objects)

    :return: a sequence of aligned SAMRecords

    """
    unaligned, widowed, paired, total = 0, 0, 0, 0
    for pair in paired_reads:
        r1_hash, r2_hash = hash_transform(index, pair)
        r1_aln, r2_aln = best_alignment(ref, index, r1_hash, r2_hash, max_edit_dist,
                                        max_indels, min_seeds, max_hits)
        if r1_aln and r2_aln:
            paired += 1
        if r1_aln or r2_aln:
            widowed += 1
        else:
            unaligned += 1
        total += 1
        try:
            yield alignment_info_to_sam(pair[0], r1_aln, pair[1].id, r2_aln, read_group,
                                        is_first=True)
            yield alignment_info_to_sam(pair[1], r2_aln, pair[0].id, r1_aln, read_group,
                                        is_first=False)
        except ValueError:
            continue
        if total % 1000 == 0:
            print('aligned: {} pairs'.format(total))
    print('Aligned: {}   Widowed:  {},  Unaligned: {}'.format(paired, widowed, unaligned))


@debug
def alignment_info_to_sam(seqrecord, aln_info, mate_id, mate_info, read_group, is_first):
    """\
    Convert the internal alignment structure into an official SAM record. The reason to
    go immediately to SAM is for compatibility with other tools, and the ablity to
    write out a file that can be easily visualized. A small annoyance is that
    downstream we will have to re-infer the location of mismatch/ins/del elements
    from the cigar string; but small price to pay.

    :param seqrecord: the read (a Bio.SeqRecord.SeqRecord object)
    :param aln_info: the alignment information
    :param mate_id: the mate for the read
    :param read_group: the read group for this alignment

    :returns: a SAMRecord for the alignment

    """
    samrecord = AlignedRead()
    samrecord.qname = seqrecord.id.rsplit(':', 1)[0]
    samrecord.seq = str(seqrecord.seq).upper()
    samrecord.is_unmapped = aln_info == None
    if aln_info:
        samrecord.mapq = 255  # TODO alignment quality?
        samrecord.pos = aln_info.offset
        samrecord.tags += [("NM", aln_info.mismatches), ("RG", read_group)]
        #samrecord.cigar = [(0, len(str(seqrecord.seq)))]  # TODO allow indels at some point
        #samrecord.cigarstring = '{}M'.format(len(str(seqrecord.seq)))
        samrecord.cigarstring = aln_info.cigar
        samrecord.rname, samrecord.tid = 0, 0  # TODO deal with multiple contigs
    else:
        samrecord.tags += [("RG", read_group)]
    if mate_info:
        samrecord.mpos = mate_info.offset
        samrecord.pnext = mate_info.offset
        samrecord.rnext = 0  # TODO deal with multiple contigs
        samrecord.mate_is_reverse = mate_info.reversed
    if aln_info and mate_info:
        # proper pair: reads are pointing at each other
        if aln_info.offset < mate_info.offset:
            samrecord.is_proper_pair = mate_info.reversed and not aln_info.reversed
        else:
            samrecord.is_proper_pair = aln_info.reversed and not mate_info.reversed
        # calculate insert
        first, second = (aln_info, mate_info) if is_first else (mate_info, aln_info)
        samrecord.isize = first.offset - second.offset if first.reversed else second.offset - first.offset
    is_reverse = aln_info is not None and aln_info.reversed
    if is_reverse:
        samrecord.seq = samrecord.seq[::-1]
    is_unmapped = aln_info == None
    mate_is_unmapped = mate_info == None
    mate_is_reverse = mate_info is not None and mate_info.reversed
    is_second = not is_first
    # TODO allow unpaired reads (the 0x1 flag)
    samrecord.flag = (0x1 | 0x2 | 0x4 * is_unmapped | 0x8 * mate_is_unmapped |
                      0x10 * is_reverse | 0x20 * mate_is_reverse |
                      0x40 * is_first | 0x80 * is_second)
    if samrecord.is_unmapped:
        samrecord.tid = -1
        samrecord.pos = -1
        samrecord.cigarstring = '*'
        samrecord.cigar = []
    else:
        if not samrecord.seq:
            # cleared by PySam (this can happen for certain cigar strings)
            samrecord.seq = str(seqrecord.seq).upper()
        try:
            if sum([x[1] for x in samrecord.cigar if x[0] != 2]) != len(samrecord.seq):
                print('ERROR AT POSITION {}'.format(samrecord.pos))
                raise ValueError('Cigar {} does not fit sequence {}'.format(samrecord.cigarstring, samrecord.seq))
        except TypeError:
            print('No seq in record: {}'.format(samrecord))
            print('WTF is goin on? \n{}'.format(seqrecord))
    samrecord.qual = ''.join([chr(q + 33) for q in seqrecord._per_letter_annotations['phred_quality']])
    return samrecord
    
