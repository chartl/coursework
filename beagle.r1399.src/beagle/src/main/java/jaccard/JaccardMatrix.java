package jaccard;

import java.lang.IllegalStateException;
import javax.print.DocFlavor;
import java.io.PrintWriter;
import java.util.Arrays;


/**
 * Created by Christopher on 7/7/2015.
 */
public class JaccardMatrix {

    private final int MAX_ALLELES, N_ROWS;
    protected final int m2Offset, m3Offset, r4Offset, maxM2, maxM3;
    private final double[][] matrix;

    private final static int[] FACTORIAL = new int[]{1, 1, 2, 6, 24, 120, 720};

    private static int smallFact(int x) {
        if ( x < 0 ) {
          return 0;
        } else if ( x < FACTORIAL.length ) {
            return FACTORIAL[x];
        } else {
            return x * smallFact(x-1);
        }
    }

    protected static int smallBinom(int n, int k) {
        if (n-k < 0) {
          return 0;
        }
        return smallFact(n)/(smallFact(k)*smallFact(n-k));
    }

    private static void assert_(boolean bool, String msg) {
      if ( ! bool ) {
        throw new IllegalStateException(String.format("Assertion is FALSE!%n%s", msg));
      }
    }

    private static void assert_(boolean bool) { assert_(bool, ""); }

    protected static int seqBefore2(int i, int j, int maxAlleles) {
        assert_(i < j);
        if ( j >= maxAlleles ) {
            return smallBinom(maxAlleles, 2);
        }
        int nSeq = 0;
        for ( int t = 0; t < i; t++ ) {
            nSeq += (maxAlleles - t - 1);
        }

        if ( j-1 > i ) {
            nSeq += ( j-1 - i );
        }

        return nSeq;
    }

    protected static int seqBefore3(int i, int j, int k, int maxAlleles) {
        assert_(i < j);
        assert_(j < k);
        if ( k >= maxAlleles ) {
            return smallBinom(maxAlleles, 3);
        }
        int nSeq = 0;
        for ( int t = 0; t < i; t++ ) {
            nSeq += smallBinom(maxAlleles - t - 1, 2);
        }

        for (int t = i+1; t < j; t++) {
            nSeq += (maxAlleles - t - 1);
        }

        if ( k-1 > j ) {
            nSeq += k- 1 - j;
        }
        return nSeq;
    }

    protected static int seqBefore4(int i, int j, int k, int l, int maxAlleles) {
        //System.out.printf("--Sb4: (%d, %d, %d, %d) max=%d%n", i, j, k, l, maxAlleles);
        assert_(i < j);
        assert_(j < k);
        assert_(k < l);
        if ( l >= maxAlleles ) {
            return smallBinom(maxAlleles, 4);
        }
        int nSeq = 0;
        for ( int t = 0; t < i; t++ ) {
            //System.out.printf("ti=%d%n", t);
            nSeq += smallBinom(maxAlleles - t - 1, 3);
        }

        for ( int t=i+1; t < j; t++) {
            //System.out.printf("tj=%d%n", t);
            nSeq += smallBinom(maxAlleles - t - 1, 2);
        }

        for ( int t = j+1; t < k; t++ ) {
            //System.out.printf("tk=%d%n", t);
            nSeq += (maxAlleles - t - 1);
        }

        if ( l - 1 > k ) {
            //System.out.printf("l-1=%d   k=%d%n", l-1, k);
            nSeq += l - 1 - k;
        }

        return nSeq;
    }

    protected static int m2Row(int i, int j, int k, int l) {
        // there are 2 alleles, so set({i, j, k, l}).size() = 2
        //System.out.printf("m2Row(%d, %d, %d, %d)%n", i, j, k, l);
        if ( i == j ) {
            if (k == l) {
                if (j < k) {
                    // xx,yy
                    return 2;
                } else {
                    // yy, xx
                    return 3;
                }
            } else {
                if ( j > k ) {
                    // yy, xy
                    return 6;
                } else {
                    // xx, xy
                    return 0;
                }
            }
        } else {
            if ( k == l ) {
                if ( j > k ) {
                    // xy, xx
                    return 1;
                } else {
                    // xy, yy
                    return 5;
                }
            } else {
                // xy, xy
                return 4;
            }
        }
    }

    protected static int m3Row(int i, int j, int k, int l) {
        //System.out.printf("m3Row(%d, %d, %d, %d)%n", i, j, k, l);
        // there are 3 alleles, so set({i, j, k, l}).size() = 3
        if ( i < k ) {
            if ( j < k ) {
                if ( k == l ) {
                    // xy zz
                    return 3;
                } else {
                    // xx yz
                    return 0;
                }
            } else {
                if ( k == l ) {
                    // xz yy
                    return 5;
                } else if ( j == k ) {
                    // xy yz
                    return 2;
                } else {
                    // xz yz
                    return 6;
                }
            }
        } else {
            if ( j > l ) {
                if ( k == l ) {
                    // yz xx
                    return 8;
                } else if (i == l ) {
                    // yz xy
                    return 9;
                } else if ( i == j ) {
                    // zz xy
                    return 11;
                } else {
                    // xz xy
                    return 4;
                }
            } else {
                if ( i == j ) {
                    // yy xz
                    return 7;
                } else if (j == l) {
                    // yz xz
                    return 10;
                } else {
                    // xy xz
                    return 1;
                }
            }
        }
    }

    protected static int countRows(int numAlleles) {
        return numAlleles + 1
              + 7 * ( 1 + smallBinom(numAlleles, 2))
                + 12 * (1 + smallBinom(numAlleles, 3)) +
                + (1 + smallBinom(numAlleles, 4)); // number of rows associated with m*
    }


    protected JaccardMatrix(int alleles) {
        this.MAX_ALLELES = alleles;
        this.N_ROWS = countRows(alleles);
        matrix = new double[N_ROWS][9];
        m2Offset = MAX_ALLELES + 1;
        maxM2 = (1 + smallBinom(MAX_ALLELES, 2));
        m3Offset = m2Offset + 7 * maxM2;
        maxM3 = (1 + smallBinom(MAX_ALLELES, 3));
        r4Offset = m3Offset + 12 * maxM3; 
        assert_(r4Offset + smallBinom(MAX_ALLELES, 4) + 1 == this.N_ROWS,
                String.format("Number of rows is bad: expected %d observed %d",
                        r4Offset + smallBinom(MAX_ALLELES, 4) + 1, this.N_ROWS));
      //System.out.printf("NROW=%d, m2=%d, m3=%d, r4=%d%n", N_ROWS, m2Offset, m3Offset, r4Offset);
    }

    public void dump(PrintWriter writer) {
        for ( int row = 0; row < matrix.length; row++ ) {
            StringBuilder bldr = new StringBuilder();
            for ( int col = 0; col < matrix[0].length; col++ ) {
                bldr.append(String.format("%.4f ", matrix[row][col]));
            }
            bldr.deleteCharAt(bldr.length()-1);
            writer.printf("%s%n", bldr.toString());
        }
    }

    protected void update(double[] hapFreqs) {
        //System.out.printf("Updating: %s%n", Arrays.toString(hapFreqs));
        int m2count = 0;
        int m3count = 0;
        int r4count = 0;
        for ( int h1 = 0; h1 < hapFreqs.length; h1++) {
            //System.out.printf("%d%n", h1);
            if (h1 < MAX_ALLELES) {
                addRow(h1, type1Row(hapFreqs[h1]));
                //System.out.printf("1a: updating %d%n", h1);
            } else {
                //System.out.printf("1b: updating %d%n", h1);
                addRow(MAX_ALLELES, type1Row(hapFreqs[h1]));
            }
            for (int h2 = h1 + 1; h2 < hapFreqs.length; h2++ ) {
                //System.out.printf("%d, %d%n", h1, h2);
                if ( h2 < MAX_ALLELES ) {
                    //System.out.printf("2a: updating %d%n", m2Offset + 7 * m2count);
                    addMatrix(m2Offset + 7 * (m2count++), type2Matrix(hapFreqs[h1], hapFreqs[h2]));
                } else {
                    //System.out.printf("2b: updating %d%n", m2Offset + 7 * (maxM2 - 1));
                    addMatrix(m2Offset + 7 * (maxM2 - 1), type2Matrix(hapFreqs[h1], hapFreqs[h2]));
                }
                for (int h3 = h2 + 1; h3 < hapFreqs.length; h3++ ) {
                    //System.out.printf("%d %d %d%n", h1, h2, h3);
                    if ( h3 < MAX_ALLELES ) {
                        //System.out.printf("3a: updating %d%n", m3Offset + 12 * m3count);
                        addMatrix(m3Offset + 12 * (m3count++), type3Matrix(hapFreqs[h1], hapFreqs[h2], hapFreqs[h3]));
                    } else {
                        //System.out.printf("3b: updating %d%n", m3Offset + 12 * (maxM3 - 1));
                        addMatrix(m3Offset + 12 * (maxM3 - 1), type3Matrix(hapFreqs[h1], hapFreqs[h2], hapFreqs[h3]));
                    }
                    for (int h4 = h3 + 1; h4 < hapFreqs.length; h4++ ) {
                        //System.out.printf("%d, %d, %d, %d%n", h1, h2, h3, h4);
                        if ( h4 < MAX_ALLELES ) {
                            //System.out.printf("4a: updating %d%n", r4Offset + r4count);
                            addRow(r4Offset + (r4count++), type4Row(hapFreqs[h1], hapFreqs[h2], hapFreqs[h3],
                                    hapFreqs[h4]));
                        } else {
                            //System.out.printf("4b: updating %d%n", r4Offset + smallBinom(MAX_ALLELES, 4));
                            addRow(r4Offset + smallBinom(MAX_ALLELES, 4), type4Row(hapFreqs[h1], hapFreqs[h2], hapFreqs[h3],
                                    hapFreqs[h4]));
                        }
                    }
                }
            }
        }
    }

    private void addRow(int idx, double[] values) {
        for (int j = 0; j < values.length; j++ ) {
            matrix[idx][j] += values[j];
        }
    }

    private void addMatrix(int firstRow, double[][] values) {
        for (int i = 0; i < values.length; i++ ) {
            addRow(firstRow + i, values[i]);
        }
    }

    private double[] type1Row(double freq) {
        double f2 = freq * freq;
        double f3 = f2 * freq;
        double f4 = f3 * freq;
        return new double[]{freq, f2, f2, f3, f2, f3, f2, f3, f4};
    }

    private double[][] type2Matrix(double p, double q) {
        double pq = p * q;
        double pq2 = pq * q;
        double p2q = pq * p;
        double p2q2 = p2q * q;
        double pq3 = pq2 * q;
        double p3q = p2q * p;
        return new double[][] {
                {0.0, 0.0, pq, 2*p2q, 0.0, 0.0, 0.0, p2q, 2 * p3q},
                {0.0, 0.0, 0.0, 0.0, pq, 2 * p2q, 0.0, p2q, 2*p3q},
                {0.0, pq, 0.0, pq2, 0.0, p2q, 0.0, 0.0, p2q2},
                {0.0, pq, 0.0, p2q, 0.0, pq2, 0.0, 0.0, p2q2},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2 * pq, pq2 + p2q, 4 * p2q2},
                {0.0, 0.0, 0.0, 0.0, pq, 2 * pq2, 0.0, pq2, 2 * pq3},
                {0.0, 0.0, pq, 2 * pq2, 0.0, 0.0, 0.0, pq2, 2 * pq3}
        };
    }

    private double[][] type3Matrix(double p, double q, double r) {
        double pqr = p * q * r;
        double p2qr = pqr * p;
        double pq2r = pqr * q;
        double pqr2 = pqr * r;
        return new double[][] {
                {0.0, 0.0, 0.0, 2*pqr, 0.0, 0.0, 0.0, 0.0, 2*p2qr},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*p2qr},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*pq2r},
                {0.0, 0.0, 0.0, 0.0, 0.0, 2 * pqr, 0.0, 0.0, 2 * pqr2},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*p2qr},
                {0.0, 0.0, 0.0, 0.0, 0.0, 2*pqr, 0.0, 0.0, 2*pq2r},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*pqr2},
                {0.0, 0.0, 0.0, 2*pqr, 0.0, 0.0, 0.0, 0.0, 2*pq2r},
                {0.0, 0.0, 0.0, 0.0, 0.0, 2*pqr, 0.0, 0.0, 2*p2qr},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*pq2r},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pqr, 4*pqr2},
                {0.0, 0.0, 0.0, 2*pqr, 0.0, 0.0, 0.0, 0.0, 2*pqr2}
        };
    }

    private double[] type4Row(double p, double q, double r, double s) {
        return new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6*p*q*r*s};
    }

    protected double[][] getMatrix() { return matrix; }
}
