"""\
Programming assignment 1

"""
from argparse import ArgumentParser
from itertools import islice
import os
from pprint import pprint

from pysam import Samfile
from vcf import Reader as VCFReader, Writer as VCFWriter

from m260b.align.hash_reference import build_hashtable, verify_hash, load_hash, save_hash
from m260b.align.hash_reads import align_pairs
from m260b.traversal.active_window import active_pileups 
from m260b.align.utils import compsam
from m260b.genotype.snp import assignment1_call_snps as call_variants 
from m260b.io.utils import read_basic_fasta, read_paired_fasta, vcf2m260
import m260b.debug.debug


SAMPLE_NAME = 'hw1_sample'


def parse_args():
    parser = ArgumentParser('M260B -- Programming Assignment 1')
    parser.add_argument('reference_file', help='Path to the assignment 1 reference')
    parser.add_argument('reads_file', help='Path to the assignment 1 reads file')
    parser.add_argument('--kmer', help='Size of the kmer to use for hashing', type=int, default=15)
    parser.add_argument('--stride', help='Stride of the hash table', type=int, default=1)
    parser.add_argument('--ref_idx', help='Path to a pre-computed reference hash table', default=None)
    parser.add_argument('--debug', help='Turn on debugging (SLOW)', action='store_true')
    parser.add_argument('--out_bam', help='Write out the intermediate bam file to this location')
    parser.add_argument('--out_vcf', help='Write out the intermediate variants to this vcf')
    parser.add_argument('--no_format', help='Omit extra format info from VCF (for IGV)', action='store_true')
    return parser.parse_args()


def get_sorted_aligned_reads(args, header, sequence):
    if args.ref_idx and os.path.exists(args.ref_idx):
        print("Loading index...")
        ref_index = load_hash(args.ref_idx)
    else:
        print("Computing reference index...")
        ref_index = build_hashtable(sequence, args.kmer, args.stride)
    print("Verifying hash...")
    for hash_, offset_ in islice(ref_index[0].iteritems(), 20):
        if not verify_hash(sequence, offset_, args.kmer, hash_):
            raise ValueError('Index failed to verify: offset {} has mismatching hashes'.format(offset_))
    print("Aligning reads...")
    pair_iterator = read_paired_fasta(args.reads_file)
    sam_iterator = align_pairs(sequence, ref_index, pair_iterator, 'hw1_rg')
    print('Sorting SAMRecords in memory...')
    sam_iterator = iter(sorted(sam_iterator, cmp=compsam))
    if args.out_bam:
         header = {'HD': {'VN': '1.0'},
                   'SQ': [{'SN': header[1:], 'LN': len(sequence)}],
                   'RG': [{'ID': 'hw1_rg', 'SM': SAMPLE_NAME, 'PU': 'Unknown', 'PL': 'Unknown', 
                          'LB': 'Unknown'}]}
         outfile = Samfile(args.out_bam, 'wb', header=header)
         for read in sam_iterator:
             outfile.write(read)
         outfile.close()
         infile = Samfile(args.out_bam, 'rb')
         sam_iterator = infile
    return sam_iterator


def main(args):
    m260b.debug.debug.DEBUG = args.debug
    ref_header, ref_sequence = read_basic_fasta(args.reference_file)
    reads = get_sorted_aligned_reads(args, ref_header, ref_sequence)
    vcf_stream = VCFWriter(open(args.out_vcf, 'wb'), make_vcf_header(args)) if args.out_vcf else None
    chr = ref_header[1:].strip()
    for variant in (call_variants(pileup, sample=SAMPLE_NAME, no_format=args.no_format)
                    for pileup in active_pileups(reads, ref_sequence, chr, flank=25, start=0, dfrac=0.5)):
        if variant is not None:
            print(vcf2m260(variant).strip())
            if vcf_stream:
                vcf_stream.write_record(variant)
    if vcf_stream:
        vcf_stream.flush()
        vcf_stream.close()
           

def make_vcf_header(args):
    """make a vcf header. PyVCF is retarded and needs a VCFReader to make a VCFWriter"""
    with open('.vcf_header.vcf', 'w') as tvcf:
        tvcf.write('##fileformat=VCFv4.1\n')
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
