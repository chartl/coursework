package jaccard;

import beagleutil.Samples;
import haplotype.BitHapPair;
import haplotype.HapPair;
import jacquard.Jacquard;
import jacquard.JacquardInference;
import jacquard.JacquardMatrix;
import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import main.NuclearFamilies;
import vcf.Marker;
import vcf.Markers;

import java.util.Arrays;
import java.util.*;

/**
 * Created by chartl on 7/8/2015.
 */

public class JaccardUnitTester {

    public static void assert_(boolean bool, String msg) {
      if ( ! bool ) {
        throw new java.lang.IllegalStateException(String.format("Assertion Failed!%n%s", msg));
      }
      System.out.print(".");
    }

    public static void assert_(boolean bool) { assert_(bool, ""); }

    public static void assertIntEq(int a, int b) {
        assert_(a == b, String.format("%d != %d", a, b));
    }

    public static void assertMatEq(double[][] expected, double[][] observed, double tol) {
        assert_(expected.length == observed.length,
                String.format("First dimension of matrices not equal (%d != %d)", expected.length, observed.length));
        assert_(expected.length > 0, "Matrix must have at least one row");
        assert_(expected[0].length == observed[0].length, "Matrices do not have the same number of columns");
        for ( int i = 0; i < expected.length; i ++ ) {
            for ( int j = 0; j < expected[i].length; j ++ ) {
                assert_(Math.abs(expected[i][j] - observed[i][j]) < tol,
                        String.format("E[%d,%d] != O[%d,%d] (%f != %f)", i, j, i, j,
                                expected[i][j], observed[i][j]));
            }
        }
    }

    private static void testUpdate() {
        JaccardMatrix matrix = new JaccardMatrix(2);
        // 2 type-A (PLUS catchall), 1 type-B (PLUS) catchall, 0 type-C  and type-d (PLUS) catchall
        assertIntEq(matrix.getMatrix().length, (2 + 1) + (1 + 1) * 7 + (0 + 1) * 12 + (0 + 1));
        matrix.update(new double[]{0.8, 0.15, 0.05});
        double[][] expected = new double[][] {
                {0.8, 0.64, 0.64, 0.512, 0.64, 0.512, 0.64, 0.512, 0.4096},
                {0.15, 0.0225, 0.0225, 0.003375, 0.0225, 0.003375, 0.0225, 0.003375, 0.00050625},
                {0.05, 0.0025, 0.0025, 0.000125, 0.0025, 0.000125, 0.0025, 0.000125, 6.25e-06},

                {0, 0, 0.12, 0.192, 0, 0, 0, 0.096, 0.1536},
                {0, 0, 0, 0, 0.12, 0.192, 0, 0.096, 0.1536},
                {0, 0.12, 0, 0.018, 0, 0.096, 0, 0, 0.0144},
                {0, 0.12, 0, 0.096, 0, 0.018, 0, 0, 0.0144},
                {0, 0, 0, 0, 0, 0, 0.24, 0.114, 0.0576},
                {0, 0, 0, 0, 0.12, 0.036, 0, 0.018, 0.0054},
                {0, 0, 0.12, 0.036, 0, 0, 0, 0.018, 0.0054},

                {0, 0, 0.0475, 0.06625, 0, 0, 0, 0.033125, 0.0515375},
                {0, 0, 0, 0, 0.0475, 0.06625, 0, 0.033125, 0.0515375},
                {0, 0.0475, 0, 0.002375, 0, 0.033125, 0, 0, 0.00165625},
                {0, 0.0475, 0, 0.033125, 0, 0.002375, 0, 0, 0.00165625},
                {0, 0, 0, 0, 0, 0, 0.095, 0.0355, 0.006625},
                {0, 0, 0, 0, 0.0475, 0.00475, 0, 0.002375, 0.0002375},
                {0, 0, 0.0475, 0.00475, 0, 0, 0, 0.002375, 0.0002375},

                {0, 0, 0, 0.012, 0, 0, 0, 0, 0.0096},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0192},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0036},
                {0, 0, 0, 0, 0, 0.012, 0, 0, 0.0006},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0192},
                {0, 0, 0, 0, 0, 0.012, 0, 0, 0.0018},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0012},
                {0, 0, 0, 0.012, 0, 0, 0, 0, 0.0018},
                {0, 0, 0, 0, 0, 0.012, 0, 0, 0.0096},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0036},
                {0, 0, 0, 0, 0, 0, 0, 0.006, 0.0012},
                {0, 0, 0, 0.012, 0, 0, 0, 0, 0.0006},

                {0, 0, 0, 0, 0, 0, 0, 0, 0}
        };

        assertMatEq(expected, matrix.getMatrix(), 1e-6);

        matrix = new JaccardMatrix(5);
        // 5 type-A (PLUS catchall),

        matrix = new JaccardMatrix(3);
        // three type-2 matrices PLUS the catch-all, and one type-3 PLUS the catch-all, zero type-4 PLUS the catch all
        assertIntEq(matrix.getMatrix().length, (3 + 1) + (3 + 1) * 7 + (1 + 1) * 12 + (0 + 1));
        matrix.update(new double[]{0.6, 0.25, 0.1, 0.05});
        expected = new double[][] {
                {0.600000, 0.360000, 0.360000, 0.216000, 0.360000, 0.216000, 0.360000, 0.216000, 0.129600},
                {0.250000, 0.062500, 0.062500, 0.015625, 0.062500, 0.015625, 0.062500, 0.015625, 0.003906},
                {0.1, 0.01, 0.01, 0.001, 0.01, 0.001, 0.01, 0.001, 0.0001},
                {0.05, 0.0025, 0.0025, 0.000125, 0.0025, 0.000125, 0.0025, 0.000125, 6.25e-06},

                {0, 0, 0.15, 0.18, 0, 0, 0, 0.09, 0.108},
                {0, 0, 0, 0, 0.15, 0.18, 0, 0.09, 0.108},
                {0, 0.15, 0, 0.0375, 0, 0.09, 0, 0, 0.0225},
                {0, 0.15, 0, 0.09, 0, 0.0375, 0, 0, 0.0225},
                {0, 0, 0, 0, 0, 0, 0.3, 0.1275, 0.09},
                {0, 0, 0, 0, 0.15, 0.075, 0, 0.0375, 0.01875},
                {0, 0, 0.15, 0.075, 0, 0, 0, 0.0375, 0.01875},

                {0, 0, 0.06, 0.072, 0, 0, 0, 0.036, 0.0432}, // 11
                {0, 0, 0, 0, 0.06, 0.072, 0, 0.036, 0.0432},
                {0, 0.06, 0, 0.006, 0, 0.036, 0, 0, 0.0036},
                {0, 0.06, 0, 0.036, 0, 0.006, 0, 0, 0.0036},
                {0, 0, 0, 0, 0, 0, 0.12, 0.042, 0.0144},
                {0, 0, 0, 0, 0.06, 0.012, 0, 0.006, 0.0012},
                {0, 0, 0.06, 0.012, 0, 0, 0, 0.006, 0.0012},


                {0, 0, 0.025, 0.0125, 0, 0, 0, 0.00625, 0.003125}, // 18
                {0, 0, 0, 0, 0.025, 0.0125, 0, 0.00625, 0.003125},
                {0, 0.025, 0, 0.0025, 0, 0.00625, 0, 0, 0.000625},
                {0, 0.025, 0, 0.00625, 0, 0.0025, 0, 0, 0.000625},
                {0, 0, 0, 0, 0, 0, 0.05, 0.00875, 0.0025},
                {0, 0, 0, 0, 0.025, 0.005, 0, 0.0025, 0.0005},
                {0, 0, 0.025, 0.005, 0, 0, 0, 0.0025, 0.0005},

                {0, 0, 0.0475, 0.04325, 0, 0, 0, 0.021625, 0.0232625}, // 25
                {0, 0, 0, 0, 0.0475, 0.04325, 0, 0.021625, 0.0232625},
                {0, 0.0475, 0, 0.002375, 0, 0.021625, 0, 0, 0.00108125},
                {0, 0.0475, 0, 0.021625, 0, 0.002375, 0, 0, 0.00108125},
                {0, 0, 0, 0, 0, 0, 0.095, 0.024, 0.004325},
                {0, 0, 0, 0, 0.0475, 0.00475, 0, 0.002375, 0.0002375}, // 30
                {0, 0, 0.0475, 0.00475, 0, 0, 0, 0.002375, 0.0002375},

                {0, 0, 0, 0.03, 0, 0, 0, 0, 0.018}, // 32
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.036},
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.015},
                {0, 0, 0, 0, 0, 0.03, 0, 0, 0.003},
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.036},
                {0, 0, 0, 0, 0, 0.03, 0, 0, 0.0075},
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.006},
                {0, 0, 0, 0.03, 0, 0, 0, 0, 0.0075},
                {0, 0, 0, 0, 0, 0.03, 0, 0, 0.018},
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.015},
                {0, 0, 0, 0, 0, 0, 0, 0.015, 0.006},
                {0, 0, 0, 0.03, 0, 0, 0, 0, 0.003},

                {0, 0, 0, 0.0235, 0, 0, 0, 0, 0.013225}, // 44
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.02645},
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.0092},
                {0, 0, 0, 0, 0, 0.0235, 0, 0, 0.001175},
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.02645},
                {0, 0, 0, 0, 0, 0.0235, 0, 0, 0.0046},
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.00235},
                {0, 0, 0, 0.0235, 0, 0, 0, 0, 0.0046},
                {0, 0, 0, 0, 0, 0.0235, 0, 0, 0.013225},
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.0092},
                {0, 0, 0, 0, 0, 0, 0, 0.01175, 0.00235},
                {0, 0, 0, 0.0235, 0, 0, 0, 0, 0.001175},

                {0, 0, 0, 0, 0, 0, 0, 0, 0.0045}
        };

        assertMatEq(expected, matrix.getMatrix(), 1e-6);
    }

    private static List<Integer> listify(int a, int b, int c, int d) {
        return Arrays.asList(a, b, c, d);
    }

    private static Map<List<Integer>, Integer> makeRowExp(List<Integer> sortedAlleles) {
        if ( sortedAlleles.size() == 1 || sortedAlleles.size() == 4 ) {
            return null;
        }
        Map<List<Integer>, Integer> stateMap = new TreeMap<>(new Comparator() {
            @Override
            public int compare(Object o1, Object o2) {
                if (! (o1 instanceof List)) {
                    throw new IllegalStateException("Only use for lists");
                }

                if (! (o2 instanceof List)) {
                    throw new IllegalStateException("Only use for lists");
                }

                List<Integer> l1 = (List<Integer>) o1;
                List<Integer> l2 = (List<Integer>) o2;

                if ( l1.size() < l2.size() ) {
                    return -1;
                }

                for (int idx = 0; idx < l1.size(); idx++ ) {
                    int c = Integer.compare(l1.get(idx), l2.get(idx));
                    if ( c != 0 ) {
                        return c;
                    }
                }

                return 0;
            }
        });
        if ( sortedAlleles.size() == 2 ) {
            int x = sortedAlleles.get(0);
            int y = sortedAlleles.get(1);
            stateMap.put(listify(x, x, x, y), 0);
            stateMap.put(listify(x, y, x, x), 1);
            stateMap.put(listify(x, x, y, y), 2);
            stateMap.put(listify(y, y, x, x), 3);
            stateMap.put(listify(x, y, x, y), 4);
            stateMap.put(listify(x, y, y, y), 5);
            stateMap.put(listify(y, y, x, y), 6);
        } else {
            int x = sortedAlleles.get(0);
            int y = sortedAlleles.get(1);
            int z = sortedAlleles.get(2);
            stateMap.put(listify(x, x, y, z), 0);
            stateMap.put(listify(x, y, x, z), 1);
            stateMap.put(listify(x, y, y, z), 2);
            stateMap.put(listify(x, y, z, z), 3);
            stateMap.put(listify(x, z, x, y), 4);
            stateMap.put(listify(x, z, y, y), 5);
            stateMap.put(listify(x, z, y, z), 6);
            stateMap.put(listify(y, y, x, z), 7);
            stateMap.put(listify(y, z, x, x), 8);
            stateMap.put(listify(y, z, x, y), 9);
            stateMap.put(listify(y, z, x, z), 10);
            stateMap.put(listify(z, z, x, y), 11);
        }

        return stateMap;
    }

    private static void testCountUpdate() {
        JaccardMatrix matrix = new JaccardMatrix(5); // this has to match HAPLO_MARGIN
        JaccardEstimator estimator = new JaccardEstimator(12);
        int n = 7;
        // generate an ordering
        int[] counts = new int[matrix.getMatrix().length];
        for ( int a = 0; a < n; a++ ) {
            for ( int b = a; b < n; b ++ ) {
                for (int c = 0; c < n; c ++ ) {
                    for ( int d = c; d < n; d ++) {
                        // (ab, cd)
                        int[] glob = new int[]{a, b, c,  d};
                        Arrays.sort(glob);
                        Integer[] w = new Integer[]{a, b, c, d};
                        HashSet<Integer> items = new HashSet<>(Arrays.asList(w));
                        List<Integer> sortedUq = new ArrayList<>(items);
                        Collections.sort(sortedUq);
                        Integer[] loc = new Integer[]{a, b, c, d};
                        List<Integer> locList = new ArrayList<>(Arrays.asList(loc));
                        int globOffset;
                        if ( items.size() == 1 ) {
                            globOffset = sortedUq.get(0) < 5 ?  sortedUq.get(0) : 5;
                        } else if (items.size() == 2) {
                            globOffset = matrix.m2Offset + JaccardMatrix.seqBefore2(sortedUq.get(0), sortedUq.get(1), 5) * 7;
                        } else if (items.size() == 3) {
                            globOffset = matrix.m3Offset + JaccardMatrix.seqBefore3(sortedUq.get(0), sortedUq.get(1), sortedUq.get(2), 5) * 12;
                        } else {
                            globOffset = matrix.r4Offset + JaccardMatrix.seqBefore4(sortedUq.get(0), sortedUq.get(1),
                                    sortedUq.get(2), sortedUq.get(3), 5);
                        }
                        Map<List<Integer>, Integer> locMap = makeRowExp(sortedUq);
                        if ( locMap != null && ! locMap.containsKey(locList)) {
                            throw new IllegalStateException(String.format("Map does not contain key %s", Arrays.toString(loc)));
                        }
                        int localOffset = locMap == null ? 0 : locMap.get(locList);
                        estimator.addGenotypeCounts(counts, glob, loc, sortedUq.size());
                        int oneIdx = 0;
                        for ( ; oneIdx < counts.length; oneIdx++) {
                            if ( counts[oneIdx] == 1 )
                                break;
                        }
                        assert_(localOffset + globOffset == oneIdx, String.format("%d + %d != %d", globOffset, localOffset, oneIdx));
                        counts[localOffset + globOffset] = 0;
                    }
                }
            }
        }
    }

    private static void testCountsBefore() {
        int mstar = 8;
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        for (int i = 0; i < 12; i++) {
            if ( i < mstar ) {
                assertIntEq(count1, i);
                count1++;
            }
            for ( int j = i + 1; j < 12; j ++) {
                if ( j < mstar ) {
                    //System.out.printf("Test SB2: (%d, %d)%n", i, j);
                    assertIntEq(JaccardMatrix.seqBefore2(i, j, mstar), count2);
                    count2++;
                } else if ( i >= mstar ) {
                    assertIntEq(JaccardMatrix.seqBefore2(i, j, mstar), count2);
                }
                for ( int k = j + 1; k < 12; k ++ ) {
                    if ( k < mstar ) {
                        //System.out.printf("Test Sb3(%d, %d, %d)%n", i, j, k);
                        assertIntEq(JaccardMatrix.seqBefore3(i, j, k, mstar), count3);
                        count3++;
                    } else if ( i >= mstar ) {
                        assertIntEq(JaccardMatrix.seqBefore3(i, j, k, mstar), count3);
                    }
                    for ( int l = k + 1; l < 12; l ++ ) {
                        if ( l < mstar ) {
                            assertIntEq(JaccardMatrix.seqBefore4(i, j, k, l, mstar), count4);
                            count4++;
                        } else if (i >= mstar) {
                            // todo: this branch is (as of now) meaningless (as are the others)
                            assertIntEq(JaccardMatrix.seqBefore4(i, j, k, l, mstar), count4);
                        }
                    }
                }
            }
        }
    }

    private static void testRowLookups() {
        int x = 1;
        int y = 2;
        int z = 3;
        int zz = 4;

        assertIntEq(JaccardMatrix.m2Row(x, x, x, y),0);
        assertIntEq(JaccardMatrix.m2Row(x, y, x, x),1);
        assertIntEq(JaccardMatrix.m2Row(x, x, y, y), 2);
        assertIntEq(JaccardMatrix.m2Row(y, y, x, x), 3);
        assertIntEq(JaccardMatrix.m2Row(x, y, x, y), 4);
        assertIntEq(JaccardMatrix.m2Row(x, y, y, y), 5);
        assertIntEq(JaccardMatrix.m2Row(y, y, x, y), 6);

        assertIntEq(JaccardMatrix.m3Row(x, x, y, z), 0);
        assertIntEq(JaccardMatrix.m3Row(x, y, x, z), 1);
        assertIntEq(JaccardMatrix.m3Row(x, y, y, z), 2);
        assertIntEq(JaccardMatrix.m3Row(x, y, z, z), 3);
        assertIntEq(JaccardMatrix.m3Row(x, z, x, y), 4);
        assertIntEq(JaccardMatrix.m3Row(x, z, y, y), 5);
        assertIntEq(JaccardMatrix.m3Row(y, y, x, z), 7);
        assertIntEq(JaccardMatrix.m3Row(x, z, y, z), 6);
        assertIntEq(JaccardMatrix.m3Row(y, z, x, x), 8);
        assertIntEq(JaccardMatrix.m3Row(y, z, x, y), 9);
        assertIntEq(JaccardMatrix.m3Row(y, z, x, z), 10);
        assertIntEq(JaccardMatrix.m3Row(z, z, x, y), 11);
    }

    private static List<HapPair> makeTrioHaps(int dad1, int dad2, int mom1, int mom2, int child1, int child2) {
        Marker mkr = new Marker(String.format("1\t400\t.\tA\tC,G,T\t100\tPASS\t.\tGT\t%d|%d\t%d|%d\t%d|%d",
                child1, child2, dad1, dad2, mom1, mom2));
        Marker[] markers = new Marker[] { mkr };
        Markers theMarkers = new Markers(markers);
        HapPair cHP = new BitHapPair(theMarkers, 1, new byte[]{(byte) child1}, new byte[]{(byte) child2});
        HapPair dHP = new BitHapPair(theMarkers, 0, new byte[]{(byte) dad1}, new byte[]{(byte) dad2});
        HapPair mHP = new BitHapPair(theMarkers, 2, new byte[]{(byte) mom1}, new byte[]{(byte) mom2});
        return Arrays.asList(new HapPair[]{dHP, cHP, mHP});
    }


    public static void testJacquardUpdate() {
        Samples samples = Samples.fromIds(new String[]{"NA12891", "NA12878", "NA12892"});
        Set<String> founders = new HashSet<>(2);
        founders.add("NA12891");
        founders.add("NA12892");
        JacquardMatrix matrix = new JacquardMatrix(samples, founders);
        matrix.update(makeTrioHaps(0, 1, 0, 1, 0, 0)); // ab aa = 4
        matrix.update(makeTrioHaps(0, 1, 0, 1, 0, 0)); // ab aa = 4
        matrix.update(makeTrioHaps(1, 0, 0, 0, 1, 0)); // ab ab = 6
        matrix.update(makeTrioHaps(0, 0, 0, 1, 0, 0)); // aa aa = 0
        matrix.update(makeTrioHaps(0, 0, 1, 0, 0, 1)); // aa ab = 2
        matrix.update(makeTrioHaps(1, 0, 1, 1, 1, 1)); // ab aa = 4
        matrix.update(makeTrioHaps(1, 0, 1, 0, 1, 1)); // ab aa = 4
        matrix.update(makeTrioHaps(1, 1, 1, 1, 1, 1)); // aa aa = 0
        matrix.update(makeTrioHaps(0, 1, 0, 1, 0, 0)); // ab aa = 4

        double[][] expected = new double[][]{
                {0, 0, 0, 0, 0.25, 0.25, 0, 0.125, 0.125},
                {0, 0, 0, 0, 0.25, 0.25, 0, 0.125, 0.125},
                {0, 0, 0, 0, 0, 0, 0.375, 0.1875, 0.140625},
                {0.75, 0.5625, 0.5625, 0.421875, 0.5625, 0.421875, 0.5625, 0.421875, 0.31640625},
                {0, 0, 0.1875, 0.28125, 0, 0, 0, 0.140625, 0.2109375},
                {0, 0, 0, 0, 0.1875, 0.28125, 0, 0.140625, 0.2109375},
                {0, 0, 0, 0, 0.25, 0.25, 0, 0.125, 0.125},
                {1, 1, 1, 1, 1, 1, 1, 1, 1},
                {0, 0, 0, 0, 0.25, 0.25, 0, 0.125, 0.125}
        };

        assert_(matrix.getPairKey(0, 1).equals("NA12891:NA12878"));

        assertMatEq(expected, matrix.getPair(0, 1), 1e-8);

        // now add in one of -every- remaining site
        matrix = new JacquardMatrix(samples, null); // no founders
        // aa bc (3)
        matrix.update(makeTrioHaps(0, 0, 2, 1, 0, 2));
        // bc aa (5)
        matrix.update(makeTrioHaps(2, 1, 0, 0, 1, 0));
        // ab ac (7)
        matrix.update(makeTrioHaps(0, 1, 0, 2, 0, 1));

        expected = new double[][]{
                {0, 0, 0, 0.05555555555555555, 0, 0, 0, 0, 0.027777777777777776},
                {0, 0, 0, 0, 0, 0.05555555555555555, 0, 0, 0.027777777777777776},
                {0, 0, 0, 0, 0, 0, 0, 0.027777777777777776, 0.05555555555555555}
        };

        assert_(matrix.getPairKey(0, 2).equals("NA12891:NA12892"));
        assertMatEq(expected, matrix.getPair(0, 2), 1e-8);
    }

    public static void testJacquardInference() {
        Matrix constrainedRows = new DenseMatrix(new double[][]{
                {0, 0.125, 0.125},
                {0, 0.125, 0.125},
                {0.375, 0.1875, 0.140625},
                {0.5625, 0.421875, 0.31640625},
                {0, 0.140625, 0.2109375},
                {0, 0.140625, 0.2109375},
                {0, 0.125, 0.125},
                {1, 1, 1},
                {0, 0.125, 0.125}
        });
        double lik = JacquardInference.logLik(constrainedRows, new DenseMatrix(new double[][]{{-0.8}, {-0.3}}));
        assert_(Math.abs(lik - 16.39926) < 1e-5,
                String.format("Likelihood doesn't match: %f != %f", lik, 16.39926));
        Matrix egrad = new DenseMatrix(new double[][]{{-1.68536},{-0.1892773}}).times(-1);
        Matrix grad = JacquardInference.grad(constrainedRows, new DenseMatrix(new double[][]{{-0.8}, {-0.3}}));
        assertMatEq(egrad.getData(), grad.getData(), 1e-3);
        double[] expected = new double[] {3.597e-6, 3.834e-3, 9.961e-1};
        double[] fit = JacquardInference.fitLikelihood(constrainedRows);
        for ( int d = 0; d < expected.length; d ++ ) {
            assert_(Math.abs(expected[d] - fit[d]) < 1e-2,
                    String.format("D%d mismatches: e(%f) != o(%f)", d, expected[d], fit[d]));
        }

        // now test the gradient empirically
        double e = 1e-3;
        double[] empJacob = new double[3];
        double[] d0 = new double[]{0.2, 0.4, 0.4};
        double baseLik = JacquardInference.simpleLik(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose());
        for ( int i = 0; i < 3; i ++ ) {
            d0[i] += e;
            double newLik = JacquardInference.simpleLik(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose());
            empJacob[i] = (newLik - baseLik)/e;
            d0[i] -= e;
        }
        double[] gradient = JacquardInference.simpleGrad(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose()).transpose().getData()[0];
        for ( int i = 0; i < gradient.length; i ++ ) {
            if ( Math.abs(empJacob[i] - gradient[i]) > 1e-2) {
                throw new IllegalStateException(String.format("Gradients do not match:%n%s%n%s",
                        Arrays.toString(empJacob),
                        Arrays.toString(gradient)));
            }
        }


        e = 1e-4;
        empJacob = new double[3];
        d0 = new double[]{0.6, 0.3, 0.1};
        baseLik = JacquardInference.simpleLik(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose());
        for ( int i = 0; i < 3; i ++ ) {
            d0[i] += e;
            double newLik = JacquardInference.simpleLik(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose());
            empJacob[i] = (newLik - baseLik)/e;
            d0[i] -= e;
        }
        gradient = JacquardInference.simpleGrad(constrainedRows, new DenseMatrix(new double[][]{d0}).transpose()).transpose().getData()[0];
        for ( int i = 0; i < gradient.length; i ++ ) {
            if ( Math.abs(empJacob[i] - gradient[i]) > 3e-3) {
                throw new IllegalStateException(String.format("Gradients do not match:%n%s%n%s",
                        Arrays.toString(empJacob),
                        Arrays.toString(gradient)));
            }
        }


        expected = new double[] {2.804e-6, 3.367e-03, 9.966e-1};
        fit = JacquardInference.fitSimplex(constrainedRows, 1e-4);
        for ( int d = 0; d < expected.length; d ++ ) {
            assert_(Math.abs(expected[d] - fit[d]) < 5e-3,
                    String.format("D%d mismatches: e(%f) != o(%f)", d, expected[d], fit[d]));
        }


    }

    public static void runTests() {
        System.out.printf("-------------- RUNNING UNIT TESTS -----------%n");
        System.out.printf("Counts before: ");
        testCountsBefore();
        System.out.printf("%nRow lookups: ");
        testRowLookups();
        System.out.printf("%nTest update: ");
        testUpdate();
        System.out.printf("%nTest counts: ");
        testCountUpdate();
        System.out.printf("%nTest Jacquard counts: ");
        testJacquardUpdate();
        System.out.printf("%nTest Jacquard inference: ");
        testJacquardInference();
        System.out.printf("All tests passed!%n");
    }
}
