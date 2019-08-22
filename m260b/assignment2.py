"""\
Programming assignment 2

Structured no differently from assignment1 except

 - Must provide a hash input/output file
 - Will use assembler rather than snp caller

"""
from argparse import ArgumentParser
from Bio.Alphabet import DNAAlphabet as DNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from itertools import islice
import os
from pprint import pprint

from pysam import Samfile
from vcf import Reader as VCFReader, Writer as VCFWriter

from m260b.align.hash_reference import build_hashtable, verify_hash, load_hash, save_hash
from m260b.align.hash_reads import align_pairs, AlignmentInfo, alignment_info_to_sam
from m260b.align.ukkonen import banded_sw 
from m260b.align.utils import compsam
from m260b.assembly.kmer_graph import build_haplotype
from m260b.calling.haplotype import vcf_from_haplotype 
import m260b.debug.debug
from m260b.io.utils import read_basic_fasta, read_paired_fasta, vcf2m260
from m260b.traversal.active_window import active_regions 


SAMPLE_NAME = 'hw2_sample'
SAM_HEADER =  lambda hdr, seq: {'HD': {'VN': '1.0'},
                              'SQ': [{'SN': hdr[1:], 'LN': len(seq)}],
                              'RG': [{'ID': 'hw2_rg', 'SM': SAMPLE_NAME, 'PU': 'Unknown', 'PL': 'Unknown', 
                                      'LB': 'Unknown'}]}

def parse_args():
    parser = ArgumentParser('M260B -- Programming Assignment 1')
    parser.add_argument('reference_file', help='Path to the assignment 1 reference')
    parser.add_argument('reads_file', help='Path to the assignment 1 reads file')
    parser.add_argument('reference_hash', help='Path to the reference hash file')
    parser.add_argument('--kmer', help='Size of the kmer to use for hashing', type=int, default=15)
    parser.add_argument('--stride', help='Stride of the hash table', type=int, default=1)
    parser.add_argument('--debug', help='Turn on debugging (SLOW)', action='store_true')
    parser.add_argument('--out_bam', help='Write out the intermediate bam file to this location')
    parser.add_argument('--haplotype_out', help='Write a synthetic bam file for the assembled haplotypes')
    parser.add_argument('--out_vcf', help='Write out the intermediate variants to this vcf')
    parser.add_argument('--no_format', help='Omit extra format info from VCF (for IGV)', action='store_true')
    parser.add_argument('--input_bam', help='Skip alignment and just run from this .bam')
    parser.add_argument('--start', help='Start at this position', type=int)
    parser.add_argument('--stop', help='Stop at this position', type=int)
    return parser.parse_args()


def get_sorted_aligned_reads(args, header, sequence):
    if args.reference_hash and os.path.exists(args.reference_hash):
        print("Loading index...")
        ref_index = load_hash(args.reference_hash)
    else:
        print("Computing reference index...")
        ref_index = build_hashtable(sequence, args.kmer, args.stride)
        save_hash(*ref_index, file=args.reference_hash)
    print("Verifying hash...")
    for hash_, offset_ in islice(ref_index[0].iteritems(), 20):
        if not verify_hash(sequence, offset_, args.kmer, hash_):
            raise ValueError('Index failed to verify: offset {} has mismatching hashes'.format(offset_))
    print("Aligning reads...")
    pair_iterator = read_paired_fasta(args.reads_file)
    sam_iterator = align_pairs(sequence, ref_index, pair_iterator, 'hw2_rg')
    sam_iterator = iter(sorted(sam_iterator, cmp=compsam))
    if args.out_bam:
         outfile = Samfile(args.out_bam, 'wb', header=SAM_HEADER(header, sequence))
         for read in sam_iterator:
             outfile.write(read)
         outfile.close()
         infile = Samfile(args.out_bam, 'rb')
         sam_iterator = infile
    return sam_iterator


def main(args):
    m260b.debug.debug.DEBUG = args.debug
    ref_header, ref_sequence = read_basic_fasta(args.reference_file)
    if args.input_bam:
        reads = Samfile(args.input_bam)
        if args.start and args.stop:
            reads = reads.fetch(ref_header[1:].strip(), args.start, args.stop)
    else:
        reads = get_sorted_aligned_reads(args, ref_header, ref_sequence)
    #vcf_stream = VCFWriter(open(args.out_vcf, 'wb'), make_vcf_header(args)) if args.out_vcf else None
    chr = ref_header[1:].strip()
    fail_reasons = Counter()
    haplo_out = None
    if args.haplotype_out:
        haplo_out = Samfile(args.haplotype_out, 'wb', header=SAM_HEADER(ref_header, ref_sequence))
    vcf_stream = VCFWriter(open(args.out_vcf, 'wb'), make_vcf_header(args)) if args.out_vcf else None 
    for region, reads in active_regions(reads, ref_sequence, chr, start_offset=0, flank=30, dfrac=1.0):
        #print('Calling region {}-{}'.format(region.start, region.stop))
        haplotype = build_haplotype(region.reference, reads, k=11, min_kmer_count=2)
        if haplotype.fail_reason:
            print('Failure {} at window\n{}'.format(haplotype.fail_reason, region))
            continue
        # align the haplotype to the reference sequence
        offset, cigar, score, mismatch = banded_sw(region.reference, haplotype.seq)
        haplotype_start = region.start + offset
        _info = AlignmentInfo(haplotype_start, cigar, False, mismatch)
        haplo_seq = SeqRecord(Seq(haplotype.seq, DNA), id='Haplotype{}'.format(region.start))
        dict.__setitem__(haplo_seq._per_letter_annotations, 'phred_quality', [40] * len(haplotype.seq))
        haplo_read = alignment_info_to_sam(haplo_seq, _info, 'nomate', None, 'hw2_rg', False)
        if haplo_out:
            haplo_out.write(haplo_read)
        #print(haplotype)
        for variant in vcf_from_haplotype(region, haplotype, SAMPLE_NAME, chr):
            if vcf_stream:
                vcf_stream.write_record(variant)
            print(vcf2m260(variant))
    if vcf_stream:
        vcf_stream.flush()
        vcf_stream.close()
         

def make_vcf_header(args):
    """make a vcf header. PyVCF is retarded and needs a VCFReader to make a VCFWriter"""
    if not os.path.exists('.vcf_header.vcf'):
        with open('.vcf_header.vcf', 'w') as tvcf:
            tvcf.write('##fileformat=VCFv4.2\n')
            tvcf.write('##CommandLine=<ID="assignment1.py",CommandLineOptions="{}">\n'.format(str(args)))
            if not args.no_format:
                tvcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
            tvcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            if not args.no_format:
                tvcf.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
                tvcf.write('##FORMAT=<ID=PL,Number=2,Type=Integer,Description="Phred-scaled Likelihoods">\n')
            tvcf.write('##reference={}\n'.format(args.reference_file))
            tvcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(SAMPLE_NAME))
    return VCFReader(filename='.vcf_header.vcf')
    


if __name__ == '__main__':
    args = parse_args()
    main(args)
