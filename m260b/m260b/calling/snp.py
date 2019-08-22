"""\
Utilities for calling SNPs

"""
from math import log

from vcf.model import _Record as VCFRecord, _Call as VCFCall, _Substitution as SNP, \
                        make_calldata_tuple as mkformat

FORMAT = mkformat(['GT', 'GQ', 'PL', 'DP'])
GT_FORMAT = mkformat(['GT'])

def assignment1_call_pileup(pileup, qual_threshold=20):
    """\
    A very basic SNP caller from a pileup for use with assignment 1.

    Assignment 1 has reads generated from the reference (+ some SNPs or indels)
    with a 2% error rate.

    The simulation sample is haploid.

    This gives P[obs|expected] = 0.98 * I(obs==expected) + 0.02 * I(obs != expected)

    with independent draws

    :param pileup: the pileup

    :return: the base of the alternate allele, and the likelihoods and qualities
             (or None * 4 for no SNP) 

    """
    if len(pileup) == 0:
        return None, None, None, None
    if pileup.counts[pileup.reference_base] == len(pileup) - pileup.counts.get(None, 0):
        return None, None, None, None
    alleles_in_order = sorted(pileup.counts.keys(), key=lambda k: pileup.counts[k])
    best_nonref = alleles_in_order[0] if alleles_in_order[0] != pileup.reference_base \
                    else alleles_in_order[1]
    qual_ref = sum((count*log(0.98 if base == pileup.reference_base else 0.02)
                    for base, count in pileup.counts.iteritems()))
    qual_alt = sum((count*log(0.98 if base == best_nonref else 0.02)
                    for base, count in pileup.counts.iteritems()))
    if qual_alt - qual_ref > qual_threshold/10.0:
        return best_nonref, int(-10 *qual_ref), int(-10 * qual_alt), \
                 int(10*(qual_alt - qual_ref))
    return None, None, None, None


def assignment1_call_snps(pileup, qual_threshold=40, sample='assignment1_sample', no_format=False):
    """\
    Call the pileup (see `assignment1_call_pileup`) and place the resulting
    variant (if any) into a PyVCF record object 


    This is a single-sample call

    """
    alt, rq, aq, gq = assignment1_call_pileup(pileup, qual_threshold)
    if alt is None:
        return None
    chrom = pileup.contig
    min_qual = min(aq, rq)
    record = VCFRecord(chrom, 1 + pileup.position, '.', pileup.reference_base,
                       [SNP(alt)], gq, [], dict(), 'GT' if no_format else 'GT:PL:GQ:DP', [])
    if no_format:
        data = GT_FORMAT(GT='1')
    else:
        data = FORMAT(GT='1', PL=[rq-min_qual, aq-min_qual], GQ=min(gq, 99), DP=len(pileup))
    call = VCFCall(record, sample, data)
    record._sample_indexes = {sample: 0}
    record.samples = [call]
    assert len(record.samples) > 0
    return record
