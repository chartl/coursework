package jacquard;

import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import main.Parameters;
import ml.optimization.LBFGS;
import ml.optimization.LBFGSOnSimplex;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by Christopher on 7/13/2015.
 */
public class JacquardInference {

    public static Map<String, double[]> fitCoefficients(JacquardMatrix matrix, Parameters params) {
        Map<String, double[]> coefs = new HashMap<>(matrix.getNSamples() * (matrix.getNSamples() + 1)/ 2);

        for ( int sample1 = 0; sample1 < matrix.getNSamples(); sample1 ++ ) {
            for ( int sample2 = 1 + sample1; sample2 < matrix.getNSamples(); sample2 ++) {
                System.out.printf("Fitting %s...%n", matrix.getPairKey(sample1, sample2));
                Matrix Q = new DenseMatrix(matrix.getPair(sample1, sample2));
                if ( params.constrain() ) {
                    Q = Q.getSubMatrix(0, Q.getRowDimension()-1, 6, 8);
                }
                double[] fit = fitSimplex(Q, (double) params.fittol());
                if ( params.constrain() ) {
                    fit = new double[]{0, 0, 0, 0, 0, 0, fit[0], fit[1], fit[2]};
                }
                coefs.put(matrix.getPairKey(sample1, sample2), fit);
            }
        }

        return coefs;
    }

    private static Matrix logit2norm(Matrix logit) {
        double[][] norm = new double[logit.getRowDimension() + 1][1];
        double sum = 0.;
        for ( int r = 0; r < logit.getRowDimension(); r ++ ) {
            norm[r][0] = 1./(1 + Math.exp(-logit.getEntry(r, 0)));
            sum += norm[r][0];
        }

        norm[norm.length - 1][0] = 1. - sum;
        return new DenseMatrix(norm);
    }

    public static double logLik(Matrix Q, Matrix logit) {
        ////System.out.printf("D=%s%n", Arrays.toString(logit2norm(logit).transpose().getData()[0]));
        Matrix probs = Q.mtimes(logit2norm(logit)).transpose();
        ////System.out.println(Arrays.toString(probs.getData()[0]));
        double sum = 0.0;
        if ( probs.getColumnDimension() <= 1 && probs.getRowDimension() > 1) {
            throw new IllegalStateException(String.format("Probabilities have the wrong shape: (%d, %d)", probs.getRowDimension(), probs.getColumnDimension()));
        }
        for (int col = 0; col < probs.getColumnDimension(); col++ ) {
            sum -= Math.log(probs.getEntry(0, col));
        }

        return sum;
    }

    public static Matrix grad(Matrix Q, Matrix theta) {
        Matrix d = logit2norm(theta);
        Matrix probs = Q.mtimes(d).transpose(); // [(K x 9) x (9 x 1)]T = 1 x K
        // the likelihood is sum(log(probs)); this means that
        // dl/d(d_j) = dl/dp_i * dp_i/d(d_j) = 1/p_i * Q_(ij)
        // now p_i is (1 x K) and Q is (K x 9), so clearly (1/p_i)Q ~ 1 x 9 is the gradient
        double[][] invProbs = new double[probs.getRowDimension()][probs.getColumnDimension()];
        for ( int i = 0; i < invProbs.length; i ++ ) {
            for (int j =0; j < invProbs[0].length; j ++) {
                invProbs[i][j] = 1./probs.getEntry(i, j);
            }
        }
        Matrix grad = (new DenseMatrix(invProbs)).mtimes(Q);
        //System.out.println(Arrays.toString(grad.getData()[0]));

        // of course, D9 = 1 - D1 - ... - D8. Therefore
        for ( int c = 0; c < theta.getRowDimension(); c ++ ) {
            grad.setEntry(0, c, grad.getEntry(0, c) - grad.getEntry(0, theta.getRowDimension()));
        }

        ////System.out.printf("grad=(%d,%d), d=(%d,%d)%n", grad.getRowDimension(), grad.getColumnDimension(), d.getRowDimension(), d.getColumnDimension());
        grad = grad.getSubMatrix(0, 0, 0, theta.getRowDimension() - 1);
        //System.out.println(Arrays.toString(grad.getData()[0]));
        d = d.transpose().getSubMatrix(0, 0, 0, theta.getRowDimension() - 1);
        return (grad.times(d).times((d.minus(1).times(-1)))).transpose().times(-1);  // G * D * (1-D)  (logistic)

    }

    public static double simpleLik(Matrix Q, Matrix delta) {
        Matrix probs = Q.mtimes(delta);
        double sum = 0;
        for ( int r = 0; r < probs.getRowDimension(); r ++ ) {
            sum += -Math.log(probs.getEntry(r, 0));
        }

        return sum;
    }

    public static Matrix simpleGrad(Matrix Q, Matrix delta) {
        Matrix probs = Q.mtimes(delta);
        double[][] invProbs = new double[probs.getRowDimension()][probs.getColumnDimension()];
        for ( int i = 0; i < invProbs.length; i ++ ) {
            for (int j =0; j < invProbs[0].length; j ++) {
                invProbs[i][j] = 1./probs.getEntry(i, j);
            }
        }
        //System.out.printf("probs: (%d, %d), Q: (%d, %d)%n", probs.getRowDimension(), probs.getColumnDimension(),
        //        Q.getRowDimension(), Q.getColumnDimension());
        Matrix grad = (new DenseMatrix(invProbs).transpose()).mtimes(Q);
        return grad.transpose().times(-1);
    }

    public static double[] fitSimplex(Matrix Q, double tol) {
        // no bells and whistles here
        Matrix delta;
        if ( Q.getColumnDimension() == 9 ) {
            delta = new DenseMatrix(new double[][]{{0.02}, {0.02}, {0.02}, {0.02}, {0.02}, {0.02}, {0.02},
                    {0.02}, {0.84}});
        } else if ( Q.getColumnDimension() == 3 ) {
            delta = new DenseMatrix(new double[][]{{0.02}, {0.02}, {0.96}});
        } else {
            throw new IllegalStateException("Q must have column dimension 9 or 3");
        }
        double initLik = simpleLik(Q, delta);
        Matrix grad = simpleGrad(Q, delta);
        double lik = initLik;
        boolean[] flags = null;
        for ( int i = 0; i < 30000; i ++ ) {
            //System.out.printf("grad: (%d, %d), delta: (%d, %d), i: %d%n", grad.getRowDimension(),
            //        grad.getColumnDimension(), delta.getRowDimension(), delta.getColumnDimension(), i);
            flags = LBFGSOnSimplex.run(grad, lik, tol, delta);
            if ( flags[0] ) {
                break;
            }
            lik = simpleLik(Q, delta);

            if ( flags[1] ) {
                grad = simpleGrad(Q, delta);
            }

            if ( i % 1000 == 0 ) {
                System.out.printf("iter=%d, initLik=%f, lik=%f, aggregateChange=%f%n", i, initLik, lik, lik - initLik);
            }
        }

        System.out.printf("Converged with initial log-likelihood %f and final %f%n", initLik, simpleLik(Q, delta));

        return delta.transpose().getData()[0];
    }

    public static double[] fitLikelihood(Matrix Q) {
        Matrix theta;
        if ( Q.getColumnDimension() == 9 ) {
            theta = new DenseMatrix(new double[][]{{-5, -5, -5, -5, -5, -5, -5, -3}}).transpose();
        } else if ( Q.getColumnDimension() == 3 ) {
            // restricted to no inbreeding
            theta = new DenseMatrix(new double[][]{{-5, -3}}).transpose();
        } else {
            throw new IllegalArgumentException("There are 9 Jacquard indeces, and therefore should be 9 columns or 3 if assuming no inbreeding");
        }
        double lik = logLik(Q, theta);
        double initLik = lik;
        Matrix gradient = grad(Q, theta);
        //System.out.printf("lik = %f, grad = %s%n", lik, Arrays.toString(gradient.transpose().getData()[0]));
        double step = 1.;
        double stepFactor = 1./Math.log(1 + Q.getRowDimension());
        int conv = 0;
        for (int i = 0; i < 50000; i ++ ) {
            step = 1.;
            //System.out.printf("grad: (%d, %d); theta: (%d, %d)%n", gradient.getRowDimension(), gradient.getColumnDimension(),
            //        theta.getRowDimension(), theta.getColumnDimension());
            Matrix thetaNew = null;
            Matrix delNew = null;
            for (int k= 0; k < 100; k ++ ) {
                thetaNew = theta.minus(gradient.times(step));
                delNew = logit2norm(thetaNew);
                if ( delNew.getEntry(delNew.getRowDimension()-1,0) < 0 ) {
                    step *= stepFactor;
                } else {
                    break;
                }
            }
            if ( delNew.getEntry(delNew.getRowDimension()-1, 0) < 0 ) {
                throw new IllegalStateException("Error: gradient attempted to push values off the simplex");
            }

            double dif = 0.;
            Matrix del = logit2norm(theta);
            for (int c = 0; c < thetaNew.getColumnDimension() + 1; c ++ ) {
                double d = Math.abs(del.getEntry(c, 0) - delNew.getEntry(c, 0))/step;
                if ( d > dif ) {
                    dif = d;
                }
            }

            theta = thetaNew;
            lik = logLik(Q, theta);
            if ( dif < 1e-8 ) {
                // tolerance
                conv = 1;
                break;
            }
            //flags = LBFGS.run(gradient, lik, 1e-5, theta);
            //if ( flags[0] ) {
            //    break;
            //}

            gradient = grad(Q, theta);
            // if (flags[1]) {
            //    gradient = grad(Q, theta);
            // }
        }

        System.out.printf("Converged (%d) at initNegLogLik=%f, finalNegLogLik=%f%n", conv, initLik, lik);

        return logit2norm(theta).transpose().getData()[0];
    }
}
