from argparse import ArgumentParser
from collections import namedtuple, OrderedDict

import numpy as np

STD_VCF_HDR = ['##fileformat=VCFv4.1', '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
VCF_CHR_HDR = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

Marker = namedtuple('Marker', ('chr', 'pos', 'ref', 'alt', 'mutated'))


def get_args():
    parser = ArgumentParser('SimHaplo')
    parser.add_argument('vcf', help='The input ancestral haplotype (vcf)')
    parser.add_argument('out_dir', help='The output directory')
    parser.add_argument('--mu_m', help='The mutation rate (events/bp)', type=float, default=1e-8)
    parser.add_argument('--mu_r', help='The recombination rate (events/bp)', type=float, default=1e-8)
    parser.add_argument('--Ne', help='The effective population size for simulation', type=int, default=15000)
    parser.add_argument('--min_r', help='The minimum number of (distinguishable) recombinations', type=int, default=1)
    parser.add_argument('--min_rpass', help='The minimum number of passes after recombination', type=int, default=0)
    parser.add_argument('--min_pass', help='The minimum number of passes for simulation', type=int, default=10000)
    parser.add_argument('--n_samples', help='The number of samples to take from the final population', default=500)
    parser.add_argument('--seed', help='Random seed', default=31245, type=int)

    return parser.parse_args()



def main():
    args = get_args()
    np.random.seed(args.seed)
    markers, hap_size, hap_a = read_ancestral(args.vcf)
    hap_counts, hap_graph = run_sim(hap_a, hap_size, markers, args)
    sample_haps = sample_haplotypes(hap_counts, args.n_samples)
    write_samples(sample_haps, markers, args.out_dir)
    write_graph(hap_graph, hap_counts, args.out_dir)


def read_ancestral(vfile):
    start, end, markers, haplotype = -1, -1, list(), list()
    with open(vfile) as infile:
        for line in infile:
            if line[0] != '#':
                fields = line.strip().split('\t')
                markers.append(Marker(fields[0], int(fields[1]), fields[3], fields[4], False))
                if '1|1' in fields[9]:
                    haplotype.append(1)
                elif '0|0' in fields[9]:
                    haplotype.append(0)
                else:
                    raise ValueError('Invalid homozygous genotype: {}'.format(fields[9]))
                if start < 0:
                    start = markers[-1].pos
                end = markers[-1].pos
    return markers, end - start + 1, tuple(haplotype)


def hapsample(counts, sumtotal=None):
    sumtotal = sumtotal or sum(counts.values())
    idx = np.random.choice(sumtotal, 1)[0]
    total = 0
    for hap, count in counts.iteritems():
        total += count
        if total > idx:
            return hap
    return counts.keys()[-1]


def run_sim(init_hap, hap_size_bp, markers, args):
    hap_counts = OrderedDict()
    hap_counts[init_hap] = args.Ne
    passes, recomb_passes, n_recombs = 0, 0, 0
    p_mut = 1 - np.exp(hap_size_bp * np.log(1 - args.mu_m))
    p_rec = 1 - np.exp(hap_size_bp * np.log(1 - args.mu_r))
    p_event = p_mut + p_rec
    print 'mutation proability is: {}'.format(p_event)
    p_event_pass = 1 - np.exp(args.Ne * np.log(1 - p_event))
    p_cond_mut = p_mut/p_event
    haplo_graph = {0: (0, hap_counts, init_hap, [])}
    while passes < args.min_pass or recomb_passes < args.min_rpass or n_recombs < args.min_r:
        print list(enumerate(hap_counts.values()))
        # is there a mutation or recombination in this passage?
        is_event = np.random.random() < p_event_pass or args.Ne in hap_counts.values()
        if is_event:
            print 'mut at pass {}'.format(passes)
            event_idx, new_haplo, parents, is_recomb = gen_event(p_cond_mut, hap_counts, markers, args)
            if parents:
                # log to the graph
                while new_haplo in hap_counts:  # redundant
                     event_idx, new_haplo, parents, is_recomb = gen_event(0.0, hap_counts, markers, args)
                haplo_graph[len(haplo_graph)] = (passes, hap_counts, new_haplo, parents)
                if is_recomb:
                    n_recombs += 1
                hap_counts = birthdeath(hap_counts, event_idx, args.Ne)
                to_die = hapsample(hap_counts, args.Ne)
                hap_counts[to_die] -= 1
                hap_counts[new_haplo] = 1
                hap_counts = birthdeath(hap_counts, args.Ne - event_idx, args.Ne)
            else:
                hap_counts = birthdeath(hap_counts, args.Ne, args.Ne)
        else:
            hap_counts = birthdeath(hap_counts, args.Ne, args.Ne)
        passes += 1
        if n_recombs > 0:
            recomb_passes += 1

    return hap_counts, haplo_graph


def gen_event(p_mut, hap_counts, markers, args):
    # which point in the passage does it occur
    event_idx = np.random.choice(args.Ne, 1)[0]
    # is it a mutation?
    is_mut = np.random.random() < p_mut
    if not is_mut:
        parent_haplo = [hapsample(hap_counts, args.Ne), hapsample(hap_counts, args.Ne)]
        if parent_haplo[0] == parent_haplo[1]:
            # not a distinguishable recombination
            return event_idx, parent_haplo[0], None, False
        # select a crossover index -- this will result in recombination rates correlating with marker density; but w/e
        recomb_idx = np.random.choice(len(parent_haplo[0])-1, 1)[0]+1
        new_haplo = tuple([parent_haplo[0][i] for i in xrange(recomb_idx)] + 
                          [parent_haplo[1][i] for i in xrange(recomb_idx, len(parent_haplo[1]))])
        assert len(new_haplo) == len(parent_haplo[0])
        return event_idx, new_haplo, parent_haplo, True
    # mutation
    mut_idx = np.random.choice(len(markers),1)
    while markers[mut_idx].mutated:
        mut_idx = np.random.choice(len(markers),1)
    markers[mut_idx] = Marker(*(list(markers[mut_idx][:-1]) + [True]))
    parent_haplo = hapsample(hap_counts, args.Ne)
    child_haplo = tuple([x for i, x in enumerate(parent_haplo) if i < mut_idx] + 
                        [1 - parent_haplo[mut_idx]] + 
                        [x for i, x in enumerate(parent_haplo) if i > mut_idx])
    assert len(child_haplo) == len(parent_haplo)
    return event_idx, child_haplo, [parent_haplo], False
    

def fsample(cdict, total=None):
    total = total or sum(cdict.values())
    freqs = [float(v)/total for v in cdict.values()]
    x = np.random.random()
    f = 0.
    for idx, fq in enumerate(freqs):
        f += fq
        if x < f:
            return cdict.keys()[idx]
    if f < 1.:
        raise ValueError('impossible: {} {} {}'.format(f, freqs, repr(list(enumerate(cdict.values())))))
    return cdict.keys()[-1]


def birthdeath(counts, n, Ne):
    # straightforward trivial no-cleverness simulation
    for _ in xrange(n):
        counts[fsample(counts, total=Ne)] += 1
        counts[fsample(counts, total=Ne)] -= 1
    return counts


def sample_haplotypes(hap_counts, n_samples):
    return [(hapsample(hap_counts), hapsample(hap_counts)) for _ in xrange(n_samples)]


def write_samples(sample_haps, markers, out_dir):
    out_vcf = '{}/sample_output.vcf'.format(out_dir)
    with open(out_vcf, 'w') as out:
        out.write('\n'.join(STD_VCF_HDR))
        out.write('\t'.join(VCF_CHR_HDR + ['sample{}'.format(x) for x in xrange(len(sample_haps))]) + '\n')
        for i, marker in enumerate(markers):
            minfo = [marker.chr, str(marker.pos), '.', marker.ref, marker.alt, '.', '.', '.', 'GT']
            for h1, h2 in sample_haps:
                minfo.append('{}|{}'.format(h1[i],h2[i]))
            out.write('\t'.join(minfo) + '\n')


def write_graph(hap_graph, hap_counts, out_dir):
    haplotype_map = '{}/haplo.map.txt'.format(out_dir)
    haplotype_graph = '{}/haplo.graph.txt'.format(out_dir)
    hap2idx = {k: i for i, k in enumerate(hap_counts.keys())}
    with open(haplotype_map, 'w') as out:
        out.write('hap_number\thap_seq\tfinal_count\n')
        for hap, cnt in hap_counts.iteritems():
            out.write('{}\t{}\t{}\n'.format(hap2idx[hap], hap, cnt))
    with open(haplotype_graph, 'w') as out:
        out.write('event_number\tpasses\tcounts\tnew_hap\tparents\n')
        for idx in xrange(len(haplotype_graph)):
            cnts = ','.join(['{}:{}'.format(k, v) for k, v in haplotype_graph[idx][1].iteritems()])
            nh = hap2idx[haplotype_graph[idx][2]]
            pa = ','.join(map(str, [hap2idx[p] for p in haplotype_graph[idx][3]]))
            out.write('{}\t{}\t{}\t{}\t{}\n'.format(idx, haplotype_graph[idx][0], cnts, nh, pa))


if __name__ == '__main__':
    main()

