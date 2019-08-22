"""
Utilities for dealing with called (assembled) haplotypes

"""
from collections import namedtuple
from itertools import izip as zip
import re
from vcf.model import _Record as VCFRecord, _Call as VCFCall, _Substitution as SUB
from vcf.model import make_calldata_tuple as mkformat
from m260b.align.hash_reads import align_pairs, AlignmentInfo, alignment_info_to_sam
from m260b.align.ukkonen import banded_sw, full_sw

CalledEvent = namedtuple('CalledEvent', ('pos', 'bases', 'type', 'qual'))

CIGAR_REGEX = re.compile('(\d+)([MID])') 

GT_FORMAT = mkformat(['GT'])

def vcf_from_haplotype(region, haplotype, sample, contig):
    for call in split_haplotype(region, haplotype):
         ref = region.reference[call.pos - region.start]
         if len(call.bases) == 1:
             alt = [SUB(call.bases)]
         elif call.type == 'INS':
             alt = [SUB(call.bases)]
         else:
             alt = [SUB(call.bases[0])]
             ref = call.bases
         record = VCFRecord(contig, 1 + call.pos, '.', 
                            ref, alt, '.' if call.qual <= 0 else str(call.qual),
                            [], dict(), 'GT', [])
         call = VCFCall(record, sample, GT_FORMAT(GT='1'))
         record._sample_indexes = {sample: 0}
         record.samples = [call] 
         yield record


def _get_score(score, offset, event_scores, indel_adj=0):
    """\
    The events are only the first break in a path, so MNPs will look like

    [400, -1, -1]

    in this case, check if the score is -1; if so, just use the last one

    """
    s = event_scores[offset + indel_adj]
    return s if s != -1 else score


def split_haplotype(region, haplotype, max_edit_per_10bp=1):
    offset, cigar, score, mismatch = banded_sw(region.reference, haplotype.seq, not_in_ref_penalty=40)
    if offset > 0:
        offset, cigar, score, mismatch = full_sw(region.reference, haplotype.seq, lenient=True)
    #print('splitting: {} // mismatch={} // pos={}'.format(cigar, mismatch, region.start))
    cigar_elements = [ (int(size), oper) for size, oper in CIGAR_REGEX.findall(cigar) ]
    if (len(cigar_elements) + mismatch) / (len(haplotype.seq)/10.) > max_edit_per_10bp:
        # bad haplotype; kill the events
        print('Bad haplotype alignment at {}'.format((region.start, region.stop)))
        haplotype = haplotype.__class__('bad_alignment', seq=[], event_scores=[])
    cigar_elements = iter(cigar_elements)
    cur_size, cur_oper = next(cigar_elements)
    ref_adj = 0
    score = -2
    for h_offset, h_base in enumerate(haplotype.seq):
        #print("o={}, ho={}/{}, adj={}, cop={}, csz={}, h_b={}, r_b={}".format(offset, h_offset, len(haplotype.seq), ref_adj, cur_oper, cur_size, h_base, region.reference[offset + h_offset + ref_adj]))
        if cur_oper == 'M':
            if h_base != region.reference[offset + h_offset + ref_adj]:
                score = _get_score(score, h_offset, haplotype.event_scores)
                # SNP (or MNP but we'll split those up)
                yield CalledEvent(region.start + h_offset + ref_adj, h_base, 'SNP', score)
            cur_size -= 1
            if cur_size == 0 and h_offset < region.stop:
                # there's some next, non-M cigar string, so yield it as an event
                cur_size, cur_oper = next(cigar_elements)
                if cur_oper == 'I':
                    score = _get_score(score, h_offset, haplotype.event_scores, indel_adj=1)  # insertion is next base
                    yield CalledEvent(region.start + h_offset + ref_adj, 
                                      h_base + haplotype.seq[1+h_offset:(1+h_offset+cur_size)],
                                      'INS', score)
                elif cur_oper == 'D':
                    score = _get_score(score, h_offset, haplotype.event_scores, indel_adj=1)  # deletion is next base
                    yield CalledEvent(region.start + h_offset + ref_adj,
                                      region.reference[(h_offset+ref_adj):(ref_adj+1+h_offset+cur_size)],
                                      'DEL', score)
        else:
            if cur_oper == 'I':
                ref_adj -= 1
                cur_size -= 1
            else:
                ref_adj += cur_size
                cur_size = 0
            if cur_size == 0:
                # must be an 'M' here
                cur_size, cur_oper = next(cigar_elements)
            
        
