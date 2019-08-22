package jaccard;

import la.matrix.DenseMatrix;
import la.matrix.Matrix;
import ml.optimization.LBFGS;

/**
 * Created by Christopher on 7/13/2015.
 */
public class BoundaryRegression {

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

    protected static double squaredLoss(Matrix indep, Matrix response, Matrix logitCoefs) {
        Matrix coefs = logit2norm(logitCoefs);
        //System.out.printf("Indep: (%d, %d) x Coefs (%d, %d) - Resp (%d, %d)%n",
        //       indep.getRowDimension(), indep.getColumnDimension(), coefs.getRowDimension(),
        //        coefs.getColumnDimension(), response.getRowDimension(), response.getColumnDimension());
        Matrix diff = indep.mtimes(coefs).minus(response);
        double loss = 0.0;
        for ( double[] row : diff.getData() ) {
            for ( double d : row ) {
                loss += d * d;
            }
        }
        return loss;
    }

    protected static Matrix grad(Matrix indep, Matrix response, Matrix logitCoefs) {
        Matrix coefs = logit2norm(logitCoefs);
        Matrix mdmr = (indep.mtimes(coefs)).minus(response);
        Matrix l2grad = (((mdmr.transpose()).mtimes(indep)).times(2)).transpose();
        for ( int x = 0; x < l2grad.getColumnDimension() - 1; x ++ ) {
            l2grad.setEntry(0, x, l2grad.getEntry(0, x) - l2grad.getEntry(0, l2grad.getColumnDimension()  - 1));
        }
        //System.out.printf("l2grad: (%d, %d), submat: (%d-%d, %d-%d)%n", l2grad.getRowDimension(),
        //        l2grad.getColumnDimension(), 0, 7, 0, 0);
        l2grad = l2grad.getSubMatrix(0, 7, 0, 0);
        //System.out.printf("coefs: (%d, %d), submat: (%d-%d, %d-%d)%n",
        //        coefs.getRowDimension(), coefs.getColumnDimension(), 0, 7, 0, 0);
        coefs = coefs.getSubMatrix(0, 7, 0, 0);
        Matrix oneMinusCoefs = (coefs.minus(1)).times(-1);
        Matrix cgrad = coefs.times(oneMinusCoefs);  // logit transform as a boundary function: d*(1-d)
        return l2grad.times(cgrad);
    }

    public static double[] train(double[][] indepMatrix, double[] response) {
        Matrix J = new DenseMatrix(indepMatrix);
        Matrix y = new DenseMatrix(new double[][] { response }).transpose();
        double[] initVals = new double[indepMatrix[0].length - 1];
        for ( int x = 0; x < initVals.length - 1; x ++ ) {
            initVals[x] = -2.079442; // slightly lower than 1/9
        }
        Matrix theta = new DenseMatrix(new double[][]{initVals}).transpose();

        //System.out.printf("J=%s, y=%s, theta=%s%n", J.toString(), y.toString(), theta.toString());
        double l2loss = squaredLoss(J, y, theta);
        Matrix gradient = grad(J, y, theta);
        double tol = 1e-7;
        boolean flags[] = null;
        for ( int iter = 0; iter < 1000; iter ++) {
            //System.out.printf("iter=%d, loss=%f%n", iter, l2loss);
            flags = LBFGS.run(gradient, l2loss, tol, theta);
            if ( flags[0] ) {
                break;
            }
            l2loss = squaredLoss(J, y, theta);
            if ( flags[1] ) {
                gradient = grad(J, y, theta);
            }
        }

        return logit2norm(theta).transpose().getData()[0];


    }
}
