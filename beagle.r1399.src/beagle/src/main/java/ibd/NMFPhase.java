package ibd;

import beagleutil.Samples;
import haplotype.BitHapPair;
import haplotype.HapPair;
import ml.clustering.L1NMF;
import ml.options.L1NMFOptions;
import vcf.Marker;
import vcf.Markers;
import vcf.VcfEmission;
import vcf.VcfWindow;

import java.util.ArrayList;
import java.util.List;

/**
 * Nonnegative matrix factorization phasing
 */
public class NMFPhase {

    private final int maxHaps;
    private final double markerThrehsold;
    public NMFPhase(int maxHaplos, double haploThresh) {
        maxHaps = maxHaplos;
        markerThrehsold=haploThresh;
    }

    public List<HapPair> phase(Samples samples, VcfEmission[] window) {
        double[][] genos = new double[window.length][samples.nSamples()];
        int marker = 0;
        double meanGT = 0.0;
        int gSize = window.length * samples.nSamples();
        Marker[] markers = new Marker[window.length];
        for ( VcfEmission e : window ) {
            assert e.marker().nAlleles() == 2;
            for ( int s = 0; s < samples.nSamples(); s++ ) {
                assert ! e.isPhased(s);
                genos[marker][s] = e.gl(s, (byte) 0, (byte) 1) + 2 * e.gl(s, (byte)1, (byte)1);
                meanGT += genos[marker][s];
            }
            markers[marker++] = e.marker();
        }
        meanGT /= gSize;
        // calc the variance
        double var = 0.0;
        for ( double[] row : genos ) {
            for ( double val : row ) {
                var += (val - meanGT) * (val - meanGT);
            }
        }
        var /= (gSize * (gSize - 1));
        L1NMFOptions options = new L1NMFOptions();
        // at r2 = 0.8, the residual var is 0.2 * var; so the typical residual looks likes sqrt(var/5)
        options.epsilon = 1e-3;
        options.gamma = gSize * Math.sqrt(0.2 * var) / (maxHaps * samples.nSamples());
        options.mu = gSize * Math.sqrt(0.2 * var) / (maxHaps * window.length);
        options.nClus = maxHaps;
        L1NMF nmf = new L1NMF(options);
        nmf.feedData(genos);
        nmf.clustering();
        int[][] assignments = new int[samples.nSamples()][2];
        byte[][] haplotypes = new byte[markers.length][maxHaps];
        for ( int i = 0; i < markers.length; i ++ ) {
            for ( int j = 0; j < maxHaps; j ++ ) {
                haplotypes[i][j] = (byte) (nmf.getIndicatorMatrix().getEntry(i, j) > markerThrehsold ? 1 : 0);
            }
        }
        for ( int m =  0; m < maxHaps; m++ ) {
            for ( int s = 0; s < samples.nSamples(); s ++ ) {
                if ( nmf.getCenters().getEntry(m, s) > nmf.getCenters().getEntry(assignments[s][0], s)) {
                    assignments[s][1] = assignments[s][0];
                    assignments[s][0] = m;
                } else if ( nmf.getCenters().getEntry(m, s) > nmf.getCenters().getEntry(assignments[s][1], s)) {
                    assignments[s][1] = m;
                }
            }
        }
        Markers fullmarkers = new Markers(markers);
        List<HapPair> haplogenos = new ArrayList<HapPair>(samples.nSamples());
        byte[] firstAllele = new byte[markers.length];
        byte[] secondAllele = new byte[markers.length];
        for (int s = 0; s < samples.nSamples(); s ++ ) {
            firstAllele = new byte[markers.length];
            secondAllele = new byte[markers.length];
            for ( int m = 0; m < markers.length; m ++ ) {
                firstAllele[m] = haplotypes[m][assignments[s][0]];
                secondAllele[m] = haplotypes[m][assignments[s][1]];
            }
            haplogenos.add(new BitHapPair(fullmarkers, s, firstAllele, secondAllele));
        }

        return haplogenos;
    }

}
