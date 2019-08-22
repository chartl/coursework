from argparse import ArgumentParser
from collections import defaultdict

import vcf

from m260b.io.utils import vcf2m260

def get_args():
    parser = ArgumentParser('Convert a list of VCFs to the m260b output')
    parser.add_argument('vcf', help='The VCF file (or a .list file containing one vcf file per line, in order)')
    parser.add_argument('genome', help='The genome line')
    parser.add_argument('out', help='The output file for the events')
    parser.add_argument('--filter_near_dels', help='Remove SNPs just before a deletion', action='store_true')
    parser.add_argument('--filter_bad_mnps', help='Remove 1I1M1D events (mismodeled MNPs)', action='store_true')
    return parser.parse_args()


def parse_vcf(vcf_path, filter_near_dels, filter_mnps, event_dict=None):
    print('Parsing {}'.format(vcf_path))
    event_dict = event_dict or defaultdict(list)
    last_ctype, last_var, last_call = None, None, None
    for i, record in enumerate(vcf.Reader(filename=vcf_path)):
        calltype, call = vcf2m260(record).split('\n')[:2]
        bad_snp = (last_ctype is not None) and (last_ctype == '>SNP' and calltype == '>DEL' and
                             record.POS - last_var.POS < 2)
        bad_mnp = (last_ctype is not None) and (last_ctype == '>INS' and calltype == '>DEL' and
                             record.POS - last_var.POS == 1 and len(str(record.ALT[0])) == 1)
        if last_call is not None:
            if not (bad_snp and filter_near_dels):
                if not (bad_mnp and filter_mnps):
                    event_dict[last_ctype].append(last_call)
                else:
                    calltype, call, record = None, None, None
            else:
                calltype, call, record = None, None, None
        if (i+1) % 10000 == 0:
            print('Parsed {} lines'.format(i+1))
        last_ctype, last_call, last_var = calltype, call, record 
    event_dict[calltype].append(call)  # put the last one in
    return event_dict


def parse_vcf_list(vcf_list_file, filter_near_dels, filter_bad_mnps):
    event_dict = defaultdict(list)
    for vcf_file in open(vcf_list_file):
        event_dict = parse_vcf(vcf_file.strip(), filter_near_dels, filter_bad_mnps, event_dict)
    return event_dict


if __name__ == '__main__':
    args = get_args()
    if args.vcf[-4:] == 'list':
        variants = parse_vcf_list(args.vcf, args.filter_near_dels, args.filter_bad_mnps)
    else:
        variants = parse_vcf(args.vcf, args.filter_near_dels, args.filter_bad_mnps)
    with open(args.out, 'w') as out:
        out.write('>{}\n'.format(args.genome))
        for ctype in ('>SNP','>INS', '>DEL'): 
            out.write('{}\n'.format(ctype))
            for call in variants[ctype]:
                out.write('{}\n'.format(call))

