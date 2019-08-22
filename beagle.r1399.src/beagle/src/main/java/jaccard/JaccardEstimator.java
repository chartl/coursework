package jaccard;

import beagleutil.BeagleUtils;
import beagleutil.Samples;
import blbutil.Pair;
import blbutil.Utilities;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import haplotype.HapPair;
import ibd.Haplotype;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import main.Parameters;
import ml.regression.LASSO;
import ml.regression.LinearRegression;
import ml.regression.Regression;
import net.sf.samtools.util.RuntimeEOFException;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import vcf.AllData;
import vcf.Marker;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;

/**
 * Estimates the 9 condensed Jaccard coefficients from binary data.
 *
 * For PHASED data, counts are recorded.
 *
 * For UNPHASED data with a REFERENCE, data is imputed using the reference panel,
 * and subsequently counts are recorded.
 *
 * For UNPHASED data WITHOUT a reference, haplotypes are quickly inferred
 * using nonnegative matrix factorization, and counts recorded.
 *
 * At the end of traversal, the pairwise count matrices are used to estimate
 * (via the method of moments) the 9 Jaccard coefficients between each sample pair.
 *
 *
 * The results are printed to the output file.
 */
public class JaccardEstimator {

    private static int HAPLO_MARGIN = 5; // this is m*
    private static int HAPLO_ROWS = JaccardMatrix.countRows(HAPLO_MARGIN);
    private static Comparator<Pair<String, Double>> REVSORT = new Comparator<Pair<String, Double>>() {
        @Override
        public int compare(Pair<String, Double> o1, Pair<String, Double> o2) {
            return -Double.compare(o1.second(), o2.second());
        }
    };
    private Map<Pair<Integer,Integer>, int[]> counts;
    private JaccardMatrix jMatrix;

    public JaccardEstimator(int nSamples) {
        jMatrix = new JaccardMatrix(HAPLO_MARGIN);
        counts = new HashMap<>();
        for (int i = 0; i < nSamples; i++) {
            for (int j = i+1; j < nSamples; j++ ) {
                Pair<Integer, Integer> samPair = new Pair<>(i, j);
                counts.put(samPair, new int[HAPLO_ROWS]);
            }
        }
    }

    private void assert_(boolean b) {
        if ( ! b ) {
            throw new IllegalStateException("Assertion failed");
        }
    }

    public void addCounts(Samples samples, List<HapPair> haplogenotypes, Set<String> founders) {
        assert_(samples.nSamples() == haplogenotypes.size());
        // first compute the haplotype frequencies
        /**
         * IMPORTANT: only drive by the HAPLOGENOTYPE index
         */
        Map<String, Double> hapFreqs = new HashMap<>();
        Map<String, Integer> hapOrder = new HashMap<>();
        Map<Integer, String[]> sampAlleles = new HashMap<>();
        int denom = (founders == null) ? (2 * samples.nSamples()) : (2 * founders.size());
        for ( int sample = 0; sample < samples.nSamples(); sample++ ) {
            String[] alleles = Utilities.getHaplotypes(haplogenotypes.get(sample));
            sampAlleles.put(sample, alleles);
            if ( founders == null || founders.contains(samples.id(sample))) {
                for (String allele : alleles) {
                    if (!hapFreqs.containsKey(allele)) {
                        hapFreqs.put(allele, 0.0);
                    }
                    hapFreqs.put(allele, hapFreqs.get(allele) + 1.0 / denom);
                }
            }
        }
        Pair<String, Double>[] frequencies = new Pair[hapFreqs.size()];
        int i = 0;
        for ( Map.Entry<String, Double> entry : hapFreqs.entrySet() ) {
            frequencies[i++] = new Pair<String, Double>(entry.getKey(), entry.getValue());
        }
        Arrays.sort(frequencies, REVSORT);
        double[] fPrimitive = new double[frequencies.length];
        i=0;
        for ( Pair<String, Double> elem : frequencies ) {
            hapOrder.put(elem.first(), i);
            fPrimitive[i++] = elem.second();
        }
        jMatrix.update(fPrimitive);

        int a, b, c, d;
        HashSet<Integer> alleles; 
        int[] globOrder = new int[4];
        Integer[] samOrder = new Integer[4];
        for ( int s0 = 0; s0 < samples.nSamples(); s0++) {
            String[] s0Allele = sampAlleles.get(s0);
            /*for ( String al : s0Allele ) {
                if ( ! hapFreqs.containsKey(al) ) {
                    String mloc = String.format("%s:%d", haplogenotypes.get(0).marker(0).chrom(), haplogenotypes.get(0).marker(0).pos());
                    throw new IllegalStateException(String.format("Sample %s countains non founder allele at %s:%n%s", samples.id(s0), mloc, al));
                }
            }*/
            for ( int s1 = 1 + s0; s1 < samples.nSamples(); s1++) {
                String[] s1Allele = sampAlleles.get(s1);
                /*for ( String al : s1Allele ) {
                    if ( ! hapFreqs.containsKey(al) ) {
                        String mloc = String.format("%s:%d", haplogenotypes.get(0).marker(0).chrom(), haplogenotypes.get(0).marker(0).pos());
                        throw new IllegalStateException(String.format("Sample %s countains non founder allele at %s:%n%s", samples.id(s1), mloc, al));
                    }
                }*/
                Pair<Integer, Integer> pair = new Pair<>(s0, s1);
                globOrder[0] = hapOrder.get(s0Allele[0]);
                globOrder[1] = hapOrder.get(s0Allele[1]);
                globOrder[2] = hapOrder.get(s1Allele[0]);
                globOrder[3] = hapOrder.get(s1Allele[1]);
                if (globOrder[0] <= globOrder[1]) {
                  samOrder[0] = globOrder[0];
                  samOrder[1] = globOrder[1];
                } else {
                  samOrder[1] = globOrder[0];
                  samOrder[0] = globOrder[1];
                }
                if ( globOrder[2] <= globOrder[3] ){
                  samOrder[2] = globOrder[2];
                  samOrder[3] = globOrder[3];
                } else {
                  samOrder[3] = globOrder[2];
                  samOrder[2] = globOrder[3];
                }
                Arrays.sort(globOrder);
                alleles= new HashSet<>(4);
                alleles.addAll(Arrays.asList(samOrder));
                addGenotypeCounts(counts.get(pair), globOrder, samOrder, alleles.size());
            }
        }

        hapFreqs.clear();
        sampAlleles.clear();
    }

    /**
     * Update the haplo-genotype count array, given the genotype state (ij, kl)
     *
     * @param pairCnt - the count array
     * @param globSort - the alleles sorted i < j < k < l
     * @param samSort - the alleles sorted within sample (a < c, b < d)
     * @param nAlleles - the total number of unique alleles (e.g. {i, j, k, l}.size())
     */
    protected void addGenotypeCounts(int[] pairCnt, int[] globSort, Integer[] samSort, int nAlleles) {
        //System.out.printf("%s%n", Arrays.toString(samSort));
        switch ( nAlleles ) {
            case 1:
                // type 1, the offset is just i
                int idx = samSort[0] < HAPLO_MARGIN ? samSort[0] : HAPLO_MARGIN;
                pairCnt[idx] += 1;
                return;
            case 2:
                // type 2, the number of previous type-2 matrices is
                int nmat = JaccardMatrix.seqBefore2(globSort[0], globSort[3], HAPLO_MARGIN);
                int row = JaccardMatrix.m2Row(samSort[0], samSort[1], samSort[2], samSort[3]);
                //System.out.printf("nmat=%d, row=%d, idx=%d%n", nmat, row, jMatrix.m2Offset + nmat * 7 + row);
                pairCnt[jMatrix.m2Offset + nmat * 7 + row] += 1;
                return;
            case 3:
                // type 3, the number of previous type-3 matrices is
                int a, b, c;
                a = globSort[0];
                b = globSort[1];
                if ( a == b ) {
                 b = globSort[2];
                 c = globSort[3];
                } else {
                 c = globSort[3];
                }
                nmat = JaccardMatrix.seqBefore3(a, b, c, HAPLO_MARGIN);
                row = JaccardMatrix.m3Row(samSort[0], samSort[1], samSort[2], samSort[3]);
                //System.out.printf("Glob=%s, Loc=%s, nmat=%d, row=%d%n", Arrays.toString(globSort), Arrays.toString(samSort), nmat, row);
                pairCnt[jMatrix.m3Offset + nmat * 12 + row] += 1;
                return;
            case 4:
                pairCnt[jMatrix.r4Offset + JaccardMatrix.seqBefore4(globSort[0],
                         globSort[1], globSort[2], globSort[3], HAPLO_MARGIN)] += 1;
                return;
            default:
                throw new IllegalStateException("Invalid case");
        }
    }

    public static void estimateAncestry(AllData data, Parameters params) {
        boolean hasRef = data.canAdvanceRefWindow();
        boolean hasTarget = data.canAdvanceTargetWindow();
        JaccardEstimator estimator;
        Samples samples;
        if ( hasRef && ! hasTarget ) {
            estimator = new JaccardEstimator(data.nRefSamples());
            samples = computeCounts(data, params, estimator);
        } else if ( hasRef ) {
            // phase the target
            estimator = new JaccardEstimator(data.nRefSamples() + data.nNonRefSamples());
            samples = phaseAndComputeCounts(data, params, estimator);
        } else {
            estimator = new JaccardEstimator(data.nNonRefSamples());
            samples = inferHapsAndComputeCounts(data, params, estimator);
        }

        estimateMoments(samples, estimator, params);
    }

    public static Samples phaseAndComputeCounts(AllData data, Parameters params, JaccardEstimator estimator) {
        throw new NotImplementedException();
    }

    public static Samples inferHapsAndComputeCounts(AllData data, Parameters params, JaccardEstimator estimator) {
        throw new NotImplementedException();
    }

    public static void estimateMoments(Samples samples, JaccardEstimator estimator, Parameters params) {
        PrintWriter writer;
        try {
            writer = new PrintWriter(new PrintStream(params.out() + ".jmt"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
        // write the jaccard matrix
        estimator.jMatrix.dump(writer);
        writer.close();
        try {
            writer = new PrintWriter(new PrintStream(params.out() + ".cnt"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
        // write the pairwise counts and coefficients
        PrintWriter jaccardWriter;
        try {
            jaccardWriter = new PrintWriter(new PrintStream(params.out()));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
        jaccardWriter.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                "Pair", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "I1", "I2", "K", "F"));
        int nSam = samples.nSamples();  // smaples.nSamples()
        int nPairs = nSam * (1 + nSam)/2;
        //System.out.printf("Calculating Jacquard for %d sample pairs%n", nPairs);
        int pcnt = 0;
        for ( int s1 = 0; s1 < nSam; s1++ ) { 
            for ( int s2 = s1 + 1; s2 < nSam; s2++ ) { 
                if ( pcnt % 1000 == 0 ) {
                    //System.out.printf("Done: %d pairs (%.2f%%)%n", pcnt, 100 * pcnt/((float) nPairs));
                }
                Pair<Integer, Integer> sampair = new Pair<Integer, Integer>(s1, s2);
                int[] counts = estimator.counts.get(sampair);
                writer.write(String.format("%s:%s", samples.id(s1), samples.id(s2)));
                for ( int c : counts ) {
                    writer.write(String.format("\t%d", c));
                }
                double[] jaccardProbs = calculateJaccardProbs(Doubles.toArray(Ints.asList(counts)),
                        estimator.jMatrix.getMatrix(), false);
                writer.write(String.format("%n"));
                double i1 = jaccardProbs[0] + jaccardProbs[1] + jaccardProbs[2] + jaccardProbs[3];
                double i2 = jaccardProbs[0] + jaccardProbs[1] + jaccardProbs[4] + jaccardProbs[5];
                double k = jaccardProbs[0] + 0.5 * (jaccardProbs[2] + jaccardProbs[4] + jaccardProbs[6]) +
                        0.25 * jaccardProbs[7];
                double f = jaccardProbs[0] + jaccardProbs[6];
                jaccardWriter.write(String.format("%s:%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f%n",
                        samples.id(s1), samples.id(s2), jaccardProbs[0], jaccardProbs[1], jaccardProbs[2],
                        jaccardProbs[3], jaccardProbs[4], jaccardProbs[5], jaccardProbs[6], jaccardProbs[7],
                        jaccardProbs[8], i1, i2, k, f));
                pcnt++;
            }
        }
        writer.close();
        jaccardWriter.close();
    }

    protected static double[] calculateJaccardProbs(double[] counts, double[][] jaccardMatrix, boolean constrainInbreeding) {
        return BoundaryRegression.train(jaccardMatrix, counts);
    }

    public static Samples computeCounts(AllData data, Parameters params, JaccardEstimator estimator) {
        int prevNHaps = -1;
        int winSize = 0;
        int window = 0;
        Set<String> founders = null;
        if ( params.getFounders() != null ) {
            try {
                founders = new HashSet<>(java.nio.file.Files.readAllLines(params.getFounders().toPath(), Charset.defaultCharset()));
                if ( founders.size() == 0 ) {
                    throw new IllegalStateException("Empty founders file.");
                }
                for ( String fId : founders ) {
                    if ( data.refSamples().index(fId) < 0 ) {
                        throw new IllegalStateException(String.format("Founder not found: \"%s\"", fId));
                    }
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
        Samples samples = data.refSamples();
        while ( data.canAdvanceRefWindow() ) {
            if (prevNHaps <= 0) {
                data.advanceRefWindow(0, 10);
                winSize = 10;
                window += 1;
                if ( window % 250 == 0 ) {
                  System.out.printf("Windows: %d, at: %d:%d%n", window,
                    data.refHaps().get(0).markers().markers()[0].chromIndex(),
                    data.refHaps().get(0).markers().markers()[0].pos());
                }
            } else {
                data.advanceRefWindow(winSize, 10 + winSize);  // just incorporate next 10 markers
                winSize += 10;
            }
            List<HapPair> refHaps = data.refHaps();
            int nHaps = BeagleUtils.countUniqueHaplotypes(refHaps);
            //System.out.printf("window=%d, size=%d, haps=%d%n", window, winSize, nHaps);
            if ( nHaps >= 4 ) {
                // sufficient data to break & compute
                estimator.addCounts(samples, refHaps, founders);
                prevNHaps = -1;
            } else if (nHaps > samples.nSamples()/12 && nHaps >= 3) {
                // too many haplotypes, break & compute
                estimator.addCounts(samples, refHaps, founders);
                prevNHaps = -1;
            } else {
                prevNHaps = nHaps;
            }
        }

        if (prevNHaps > 0 ) {
            estimator.addCounts(samples, data.refHaps(), founders);
        }

        return samples;
    }
}
