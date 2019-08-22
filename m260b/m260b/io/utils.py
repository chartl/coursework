"""\
Utilities for reading in/writing out files in the course formats. In particular, dealing with the peculiarities
of the input paired-end read fasta format (as opposed to say, fastQ)

"""
from Bio.Alphabet import DNAAlphabet as GENOMIC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_paired_fasta(fasta_path, default_q=17):
    """\
    Read in the paired-end read data from the fasta file into pairs of fastq records instead
    of raw strings.

    """
    def _seq_to_fq(seq, id):
        """Create a fastq record from just a sequence"""
        record = SeqRecord(Seq(seq, GENOMIC), id=id) 
        dict.__setitem__(record._per_letter_annotations, 'phred_quality',
                         [default_q] * len(seq))
        return record

    with open(fasta_path) as input:
        header = input.readline()
        if header[0] != '>':
            raise ValueError('Expect paired fasta to have a single fasta header line')
        pair_num=1
        pfx = header.strip()[1:]
        id_format = pfx + ':{}:{}'
        for line in input:
            seq1, seq2 = line.strip().split(',')
            yield _seq_to_fq(seq1, id_format.format(pair_num, 1)), _seq_to_fq(seq2, id_format.format(pair_num, 2))
            pair_num += 1


def read_basic_fasta(fasta_path):
    """\
    Just get the name and sequence for a 1-contig fasta file
    """
    with open(fasta_path) as fasta:
        header = fasta.readline().strip()
        sequence = ''.join(x.strip() for x in fasta)
    return header, sequence


def vcf2m260(record):
    """Convert a PyVCF._Record object to the m260b output format"""
    if record.is_snp:
        return ('>SNP\n'
                '{},{},{}\n').format(record.REF, record.ALT[0], record.POS - 1)
    if record.is_indel:
        if record.is_deletion:
            return ('>DEL\n'
                    '{},{}\n').format(record.REF[1:], record.POS)  # drop ref base
        # insertion
        return ('>INS\n'
                '{},{}\n').format(str(record.ALT[0])[1:], record.POS)
    raise ValueError('Variant type unrecognized: {}'.format(record))      
