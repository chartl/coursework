package jacquard;

import beagleutil.BeagleUtils;
import beagleutil.Samples;
import haplotype.HapPair;
import main.Parameters;
import net.sf.samtools.util.RuntimeEOFException;
import vcf.AllData;
import java.io.*;
import java.nio.charset.Charset;
import java.util.*;

import main.Parameters;
import vcf.AllData;

/**
 * Created by Christopher on 7/13/2015.
 */
public class Jacquard {

    public static void estimateCoefficients(AllData data, Parameters params) {
        boolean hasRef = data.canAdvanceRefWindow();
        boolean hasTarget = data.canAdvanceTargetWindow();
        if ( ! hasRef ) {
            throw new IllegalStateException("Jacquard requires phased haplotypes (for now)");
        }

        if ( hasTarget ) {
            throw new IllegalStateException("Jacquard requires phased haplotypes ONLY");
        }
        Samples allSamples = data.refSamples();
        Set<String> founders = getFounders(params);
        JacquardMatrix matrix = new JacquardMatrix(allSamples, founders);
        tabulateStates(data, params, matrix);
        Map<String, double[]> coefficients = JacquardInference.fitCoefficients(matrix, params);
        printCoefficients(coefficients, matrix, params);
        finish(params);
    }

    protected static Set<String> getFounders(Parameters params) {
        Set<String> founders = null;
        if ( params.getFounders() != null ) {
            try {
                founders = new HashSet<>(java.nio.file.Files.readAllLines(params.getFounders().toPath(), Charset.defaultCharset()));
                if ( founders.size() == 0 ) {
                    throw new IllegalStateException("Empty founders file.");
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        return founders;
    }

    protected static void tabulateStates(AllData data, Parameters params, JacquardMatrix matrix){
        while (data.canAdvanceRefWindow()) {
            data.advanceRefWindow(0, params.advanceSize());
            if (BeagleUtils.countUniqueHaplotypes(data.refHaps()) >= params.minNHaps()) {
                matrix.update(data.refHaps());
                if ( data.canAdvanceRefWindow() ) {
                    data.advanceRefWindow(0, params.skipSize());
                }
            } else {
                if ( data.canAdvanceRefWindow() ) {
                    data.advanceRefWindow(0, params.extendSize());
                } else {
                    break;
                }
            }
        }
    }

    protected static void printCoefficients(Map<String, double[]> coefs,JacquardMatrix matrix, Parameters params) {
        System.out.printf("Writing data%n");
        PrintWriter writer;
        if ( params.emitjmd() ) {
            try {
                writer = new PrintWriter(new PrintStream(params.out() + ".jmt"));
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }

            matrix.dump(writer);
            writer.close();
        }

        try {
            writer = new PrintWriter(new PrintStream(params.out()));
        } catch (FileNotFoundException e ) {
            throw new RuntimeEOFException(e);
        }

        writer.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",
                "Pair", "Constraint", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9",
                "I1", "I2", "F", "K", "R");
        String cons = params.constrain() ? "CONSTRAINED" : "UNCONSTRAINED";
        for ( Map.Entry<String, double[]> pairCoefs : coefs.entrySet() ) {
            double[] c = pairCoefs.getValue();
            double kinship = c[0] + 0.5*(c[2] + c[4] + c[6]) + 0.25 * c[7];
            double coancestry = 2 * kinship;
            double fraternity = c[0] + c[6];
            double i1 = c[0] + c[1] + c[2] + c[3];
            double i2 = c[0] + c[1] + c[4] + c[5];
            writer.printf("%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f%n",
                    pairCoefs.getKey(), cons, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], i1, i2, fraternity,
                    kinship, coancestry);
        }

        writer.close();
    }

    private static void finish(Parameters params) {
        // because of issues with the LAML library, it's useful to write a .done file
        try {
            PrintWriter writer = new PrintWriter(new PrintStream(params.out() + ".done"));
            writer.println("OK");
            writer.close();
        } catch (FileNotFoundException e) {
            throw new RuntimeEOFException(e);
        }
    }
}
