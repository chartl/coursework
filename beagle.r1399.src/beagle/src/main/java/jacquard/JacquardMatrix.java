package jacquard;

import beagleutil.BeagleUtils;
import beagleutil.Samples;
import haplotype.HapPair;

import java.io.PrintWriter;
import java.util.*;

/**
 * Created by Christopher on 7/13/2015.
 */
public class JacquardMatrix {

    private Map<String, List<double[]>> pairRows;
    private List<Integer> founders;
    private Samples allSamples;

    public JacquardMatrix(Samples samples, Set<String> founders) {
        allSamples = samples;
        if ( founders != null ) {
            this.founders = new ArrayList<Integer>();
            for ( String fId : founders ) {
                int fIdx = samples.index(fId);
                if ( fIdx < 0 ) {
                    throw new IllegalStateException(String.format("Founder %s not found", fId));
                }
                this.founders.add(fIdx);
            }
        } else {
            this.founders = new ArrayList<>(samples.nSamples());
            for ( int i = 0; i < samples.nSamples(); i ++ ) {
                this.founders.add(i);
            }
        }

        if ( this.founders == null || this.allSamples == null) {
            throw new IllegalStateException("Founders and samples must be non-null");
        }

        pairRows = new HashMap<>(samples.nSamples() * (1 + samples.nSamples())/2);
        for ( int s0 = 0; s0 < samples.nSamples(); s0++ ) {
            for ( int s1 = s0 + 1; s1 < samples.nSamples(); s1 ++) {
                String pair = getPairKey(s0, s1);
                pairRows.put(pair, new ArrayList<double[]>(20000));  // 7G for 200 samples
            }
        }
    }

    public int getNSamples() { return allSamples.nSamples(); }

    public String getPairKey(int sample1, int sample2) {
        return String.format("%s:%s", allSamples.id(sample1), allSamples.id(sample2));
    }

    public double[][] getPair(int sample1, int sample2) {
        if ( sample1 >= sample2 ) {
            throw new IllegalArgumentException("Sample 1 must be < Sample 2");
        }

        List<double[]> rows = pairRows.get(getPairKey(sample1, sample2));
        if ( rows.size() == 0 ) {
            return new double[0][0];
        }

        double[][] vals = new double[rows.size()][rows.get(0).length];
        int r = 0;
        for ( double[] row : rows ) {
            int c = 0;
            for ( double d : row ) {
                vals[r][c++] = d;
            }
            r++;
        }

        return vals;
    }

    public void update(List<HapPair> hapPairs) {
        GenotypeContainer genotypes = new GenotypeContainer(hapPairs, this.founders, this.allSamples);
        for ( int s0 = 0; s0 < this.allSamples.nSamples(); s0 ++ ) {
            for ( int s1 = s0 + 1; s1 < this.allSamples.nSamples(); s1 ++ ) {
                String id0 = allSamples.id(s0);
                String id1 = allSamples.id(s1);
                String pair = getPairKey(s0, s1);
                String[] g0 = genotypes.genotypes.get(id0);
                String[] g1 = genotypes.genotypes.get(id1);
                double[] f0 = genotypes.getFrequencies(g0);
                double[] f1 = genotypes.getFrequencies(g1);
                if ( f0 != null && f1 != null ) {
                    //System.out.println(pair);
                    // recombination: have some non-founder haplotype, just skip the window
                    pairRows.get(pair).add(jacquardRow(g0, g1, f0, f1));
                }
            }
        }
    }


    public static double[] jacquardRow(String[] g1, String[] g2, double[] f1, double[] f2) {
        //System.out.printf("g1=%s, g2=%s, f1=%s, f2=%s%n", Arrays.toString(g1), Arrays.toString(g2),
        //        Arrays.toString(f1), Arrays.toString(f2));
        if (g1[0].equals(g1[1])) {
            if (g2[0].equals(g2[1])) {
                if (g1[1].equals(g2[0])) {
                    // aaaa
                    //System.out.println("AA AA");
                    double p = f1[0];
                    double p2 = p * p;
                    double p3 = p * p * p;
                    return new double[]{p, p2, p2, p3, p2, p3, p2, p3, p * p3};
                } else {
                    // aa bb
                    //System.out.println("AA BB");
                    double p = f1[0];
                    double q = f2[0];
                    return new double[]{0, p * q, 0, p * q * q, 0, p * p * q, 0, 0, p * p * q * q};
                }
            } else {
                // aa ab or aa bc
                if (g1[1].equals(g2[0]) || g1[1].equals(g2[1])) {
                    // aa ab
                    //System.out.println("AA AB");
                    double p = f1[0];
                    double q = (g1[1].equals(g2[0])) ? f2[1] : f2[0];
                    return new double[]{0, 0, p * q, 2 * p * p * q, 0, 0, 0, p * p * q, 2 * p * p * p * q};
                } else {
                    // aa bc
                    //System.out.println("AA BC");
                    double p = f1[0];
                    double q = f2[0];
                    double r = f2[1];
                    return new double[]{0, 0, 0, 2 * p * q * r, 0, 0, 0, 0, 2 * p * p * q * r};
                }
            }
        } else {
            // (ab aa), (bc, aa), (ab, ab), (ab, ac), (ab, cd)
            if ( g2[0].equals(g2[1]) ) {
                // (ab aa), (bc, aa)
                if ( g1[0].equals(g2[0]) || g1[1].equals(g2[0]) ) {
                    //System.out.println("AB AA");
                    // ab aa
                    double p = f2[0];
                    double q = (g1[0].equals(g2[0])) ? f1[1] : f1[0];
                    //System.out.printf("p=%f, q=%f%n", p, q);
                    return new double[] {0, 0, 0, 0, p * q, 2* p * p * q, 0, p * p * q, 2 * p * p * p * q};
                } else {
                    //System.out.println("BC AA");
                    // bc aa
                    double p = f2[0];
                    double q = f1[0];
                    double r = f1[1];
                    return new double[] {0, 0, 0, 0, 0, 2 * p * q * r, 0, 0, 2 * p * p * q * r};
                }
            } else {
                // (ab, ab), (ab, ac), or (ab, cd)
                if ( g1[0].equals(g2[0]) && g1[1].equals(g2[1]) ) {
                    // ab ab [1]
                    //System.out.println("AB AB [1]");
                    double p = f1[0];
                    double q = f1[1];
                    return new double[] {0, 0, 0, 0, 0, 0, 2 * p * q, p * q * (p + q), 4 * p * p * q * q};
                } else if ( g1[1].equals(g2[0]) && g1[0].equals(g2[1]) ) {
                    //System.out.println("AB BA [2]");
                    double p = f1[0];
                    double q = f1[1];
                    return new double[] {0, 0, 0, 0, 0, 0, 2 * p * q, p * q * (p + q), 4 * p * p * q * q};
                    // ab ab [2]
                } else if (g1[0].equals(g2[0]) || g1[0].equals(g2[1]) ) {
                    // ab ac [1]
                    //System.out.println("AB AC [1]");
                    double p = f1[0];
                    double q = f1[1];
                    double r = (g1[0].equals(g2[0])) ? f2[1] : f2[0];
                    return new double[] {0, 0, 0, 0, 0, 0, 0, p * q * r, 4 * p * p * q * r};
                } else if (g1[1].equals(g2[0]) || g1[1].equals(g2[1]) ) {
                    // ab ac [2]
                    //System.out.println("AB CA [2]");
                    double p = f1[0];
                    double q = f1[1];
                    double r = (g1[1].equals(g2[0])) ? f2[1] : f2[0];
                    return new double[] {0, 0, 0, 0, 0, 0, 0, p * q * r, 4 * p * p * q * r};
                } else {
                    //System.out.println("AB CD");
                    // ab cd
                    return new double[] {0, 0, 0, 0, 0, 0, 0, 0, 4 * f1[0] * f2[0] * f1[1] * f2[1]};
                }
            }
        }
    }


    public void dump(PrintWriter writer) {
        for ( Map.Entry<String, List<double[]>> entry : this.pairRows.entrySet() ) {
            for ( double[] row : entry.getValue() ) {
                writer.printf("%s", entry.getKey());
                for ( double val : row ) {
                    writer.printf("\t%f", val);
                }
                writer.printf("%n");
            }
        }
    }


    protected class GenotypeContainer {
        public String[] haplotypes;
        public double[] frequencies;
        public Map<String, Double> frequencyMap;
        public Map<String, String[]> genotypes;

        private Comparator<Map.Entry<String, Double>> freqComparator = new Comparator<Map.Entry<String, Double>>() {
            @Override
            public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
                int dcmp = - Double.compare(o1.getValue(), o2.getValue());
                if ( dcmp != 0 ) {
                    return dcmp;
                }

                return o1.getKey().compareTo(o2.getKey());
            }
        };

        public GenotypeContainer(List<HapPair> genotypes, List<Integer> founders, Samples samples) {
            Map<String, Double> freqs = new HashMap<>(2 * founders.size());
            this.genotypes = new HashMap<>(samples.nSamples());
            for ( int i : founders ) {
                String[] locusGen = BeagleUtils.getHaplotypes(genotypes.get(i));
                for ( String haplotype : locusGen ) {
                    if ( ! freqs.containsKey(haplotype) ) {
                        freqs.put(haplotype, 0.);
                    }
                    freqs.put(haplotype, freqs.get(haplotype) + 1./(2 * founders.size()));
                }
                this.genotypes.put(samples.id(i), locusGen);
            }
            List<Map.Entry<String, Double>> freqEntries = new ArrayList<>(freqs.entrySet());
            Collections.sort(freqEntries, freqComparator);
            int idx = 0;
            haplotypes = new String[freqEntries.size()];
            frequencies = new double[freqEntries.size()];
            for ( Map.Entry<String, Double> e : freqEntries ) {
                haplotypes[idx] = e.getKey();
                frequencies[idx++] = e.getValue();
            }
            frequencyMap = freqs;

            for ( HapPair genotype : genotypes ) {
                String id = samples.id(genotype.idIndex());
                if ( ! this.genotypes.containsKey(id)) {
                    this.genotypes.put(id, BeagleUtils.getHaplotypes(genotype));
                }
            }
        }

        public double[] getFrequencies(String[] haplotypes) {
            double[] freqs = new double[haplotypes.length];
            for ( int i = 0; i < haplotypes.length; i ++ ) {
                if ( frequencyMap.containsKey(haplotypes[i]) ) {
                    freqs[i] = frequencyMap.get(haplotypes[i]);
                } else {
                    return null;
                }
            }

            return freqs;
        }
    }
}
