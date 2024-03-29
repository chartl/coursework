/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package main;

import blbutil.Const;
import blbutil.Validate;
import jaccard.JaccardUnitTester;
import pedsim.PedSimUnitTester;

import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code Parameters} represents the parameters for a Beagle analysis.
 * </p>
 * Instances of class {@code Parameters} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Parameters {

    private final String[] args;

    // data input/output parameters
    private final File gt;
    private final File gl;
    private final File gtgl;
    private final File ref;
    private final String out;
    private final File excludesamples;
    private final File excludemarkers;
    private final File ped;
    private final String chrom;
    private final float maxlr;

    // algorithm parameters
    private final int nthreads;
    private final int window;
    private final int overlap;
    private final boolean gprobs;
    private final boolean impute;
    private final boolean usephase;
    private final float singlescale;
    private final float duoscale;
    private final float trioscale;
    private final int burnin_its;
    private final int phase_its;
    private final int impute_its;
    private final long seed;
    private final boolean jaccard;

    // ibd parameters
    private final boolean ibd;
    private final float ibdlod;
    private final float ibdscale;
    private final int ibdtrim;
    private final File map;
    private final File founders;

    // expert parameters
    private final int nsamples;
    private final int buildwindow;

    // pedigree simulation parameters
    private final boolean simped;
    private final float hapfrq;

    // jacquard parameters
    private final boolean jacquard;
    private final int advancesize;
    private final int extendsize;
    private final int minnhaps;
    private final int skipsize;
    private final float fittol;
    private final boolean constrain;
    private final boolean emitjmt;

    /**
     * Constructs a new {@code Parameters} instance.
     * @param args the Beagle command line arguments.
     * @throws IllegalArgumentException if the command line arguments
     * are incorrectly specified.
     */
    public Parameters(String[] args) {

        int IMAX = Integer.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;

        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        // THIS IS WHAT YOU GET FOR NOT HAVING A TESTING FRAMEWORK GRRRR
        boolean unittest = Validate.booleanArg("unittest", argsMap, false, false);
        if ( unittest ) {
            JaccardUnitTester.runTests();
            PedSimUnitTester.runTests();
        }

        // data input/output parameters
        gt = Validate.getFile(
                Validate.stringArg("gt", argsMap, false, null, null));
        gl = Validate.getFile(
                Validate.stringArg("gl", argsMap, false, null, null));
        gtgl = Validate.getFile(
                Validate.stringArg("gtgl", argsMap, false, null, null));
        ref = Validate.getFile(
                Validate.stringArg("ref", argsMap, false, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        excludesamples = Validate.getFile(
                Validate.stringArg("excludesamples", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));
        ped = Validate.getFile(
                Validate.stringArg("ped", argsMap, false, null, null));
        chrom = Validate.stringArg("chrom", argsMap, false, null, null);
        maxlr = Validate.floatArg("maxlr", argsMap, false, 5000.0f, 1.1f, FMAX);

        // algorithm parameters
        window = Validate.intArg("window", argsMap, false, 50000, 1, IMAX);
        overlap = Validate.intArg("overlap", argsMap, false, 3000, 0, IMAX);
        gprobs = Validate.booleanArg("gprobs", argsMap, false, true);
        impute = Validate.booleanArg("impute", argsMap, false, true);
        usephase=Validate.booleanArg("usephase", argsMap, false, false);
        singlescale = Validate.floatArg("singlescale", argsMap, false, 0.8f, FMIN, FMAX);
        duoscale = Validate.floatArg("duoscale", argsMap, false, 1.0f, FMIN, FMAX);
        trioscale = Validate.floatArg("trioscale", argsMap, false, 1.0f, FMIN, FMAX);
        burnin_its = Validate.intArg("burnin-its", argsMap, false, 5, 0, IMAX);
        phase_its = Validate.intArg("phase-its", argsMap, false, 5, 0, IMAX);
        impute_its = Validate.intArg("impute-its", argsMap, false, 5, 0, IMAX);
        nthreads = Validate.intArg("nthreads", argsMap, false, 1, 1, 100000);
        seed = Validate.longArg("seed", argsMap, false, -99999, LMIN, LMAX);

        // ibd parameters
        ibd = Validate.booleanArg("ibd", argsMap, false, false);
        ibdlod = Validate.floatArg("ibdlod", argsMap, false, 3.0f, FMIN, FMAX);
        ibdscale = Validate.floatArg("ibdscale", argsMap, false, 0.0f, 0.0f, FMAX);
        ibdtrim = Validate.intArg("ibdtrim", argsMap, false, 40, 0, IMAX);
        jaccard = Validate.booleanArg("jaccard", argsMap, false, false);
        founders = Validate.getFile(Validate.stringArg("founders", argsMap, false, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, false, null, null));

        // expert parameters
        nsamples = Validate.intArg("nsamples", argsMap, false, 4, 1, IMAX);
        buildwindow = Validate.intArg("buildwindow", argsMap, false, 1200, 1, IMAX);

        // pedigree simulation parameters
        simped = Validate.booleanArg("simped", argsMap, false, false);
        hapfrq = Validate.floatArg("hapfrq", argsMap, false, 0.05f, 1e-4f, 0.2f);

        jacquard = Validate.booleanArg("jacquard", argsMap, false, false);
        advancesize = Validate.intArg("advancesize", argsMap, false, 5, 1, 1500);
        skipsize = Validate.intArg("skipsize", argsMap, false, 20, 1, 1500);
        extendsize = Validate.intArg("extendsize", argsMap, false, 1, 1, 1500);
        minnhaps = Validate.intArg("minnhaps", argsMap, false, 3, 3, 1500);
        fittol = Validate.floatArg("fittol", argsMap, false, 0.00005f, 0.000000001f, 0.1f);
        constrain = Validate.booleanArg("constrain", argsMap, false, false);
        emitjmt = Validate.booleanArg("emitjmt", argsMap, false, false);

        Validate.confirmEmptyMap(argsMap);
    }

    /**
     * Returns the Beagle command line arguments.
     * @return the Beagle command line arguments.
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a string summary of the Beagle command line arguments.
     * @return a string summary of the Beagle command line arguments.
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Command line syntax: java -jar beagle.jar [arguments]" + nl
                + nl
                + "data input/output parameters ..." + nl
                + "  gt=<VCF file: use GT field>                        (optional)" + nl
                + "  gl=<VCF file: use GL/PL field>                     (optional)" + nl
                + "  gtgl=<VCF file: use GT and GL/PL fields>           (optional)" + nl
                + "  ref=<VCF file with phased genotypes>               (optional)" + nl
                + "  out=<output file prefix>                           (required)" + nl
                + "  excludesamples=<file with 1 sample ID per line>    (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>    (optional)" + nl
                + "  ped=<linkage format pedigree file>                 (optional)" + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)" + nl
                + "  maxlr=<max GL/PL likelihood ratio>                 (default=5000)" + nl + nl

                + "algorithm parameters ..." + nl
                + "  nthreads=<number of threads>                       (default=1)" + nl
                + "  window=<markers per window>                        (default=50000)" + nl
                + "  overlap=<overlap between windows>                  (default=3000)" + nl
                + "  gprobs=<print GP field (true/false)>               (default=true)" + nl
                + "  impute=<impute ungenotyped variants (true/false)>  (default=true)" + nl
                + "  usephase=<use phase in \"gt\" or \"gtgl\" file>        (default=false)" + nl
                + "  singlescale=<model scale for singles>              (default=1.0)" + nl
                + "  duoscale=<model scale for duos>                    (default=1.0)" + nl
                + "  trioscale=<model scale trios>                      (default=1.0)" + nl
                + "  burnin-its=<number of iterations>                  (default=5)" + nl
                + "  phase-its=<number of iterations>                   (default=5)" + nl
                + "  impute-its=<number of iterations>                  (default=5)" + nl
                + "  seed=<random seed>                                 (default=-99999)" + nl + nl

                + "IBD parameters ..." + nl
                + "  ibd=<perform IBD detection (true/false)>           (default=false)" + nl
                + "  ibdlod=<min LOD score for reporting IBD>           (default=3.0)" + nl
                + "  ibdscale=<model scale factor for Refined IBD>      (default: data-dependent)" + nl
                + "  ibdtrim=<markers at each segment end>              (default=40)" + nl
                + "  jaccard=<estimate jaccard (true/false)>            (default=false)" + nl
                + "  founders=<file containing founder ids>             (optional)" + nl + nl

                + "Simulation parameters ..." + nl
                + "  simped=<simulate pedigree (true/false)>             (default=false)" + nl
                + "  hapfrq=<minimum hap frequency>                      (default=0.05)" + nl + nl

                + "Jacquard parameters ..." + nl
                + "  advancesize=<number of initial markers>             (default=5)" + nl
                + "  skipsize=<number of markers between windows>        (default=20)" + nl
                + "  extendsize=<n markers for hap extension>            (default=1)" + nl
                + "  minnhaps=<minimum number of haplotypes>             (default=3)" + nl
                + "  fittol=<fit tolerance>                              (default=5e-5)" + nl
                + "  constrain=<constrain to no inbreeding(t/f)>         (default=false)" + nl
                + "  emitjmt<emit the estimation matrix (t/f)>           (default=false)";
    }

    /**
     * Returns a sample-size-adjusted IBD scale parameter. Returns
     * {@code this.ibdscale()} if {@code this.ibdscale()!=0.0f}, and
     * returns {@code Math.max(2.0f, (float) Math.sqrt(nSamples/100.0))}
     * otherwise.
     *
     * @param nSamples the number of samples.
     * @return a sample-size-adjusted IBD scale parameter.
     * @throws IllegalArgumentException if {@code nSamples<0}.
     */
    public float adjustedIbdScale(int nSamples) {
        if (nSamples <= 0) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        if (ibdscale==0) {
            return Math.max(2.0f, (float) Math.sqrt(nSamples/100.0));
        }
        else {
            return ibdscale;
        }
    }

    // data input/output parameters

    /**
     * Returns the gt parameter or {@code null} if no gt parameter was
     * specified.
     * @return the gt parameter or {@code null} if no gt parameter was
     * specified.
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the gl parameter or {@code null} if no gl parameter was
     * specified.
     * @return the gl parameter or {@code null} if no gl parameter was
     * specified.
     */
    public File gl() {
        return gl;
    }

    /**
     * Returns the gtgl parameter or {@code null} if no gtgl parameter was
     * specified.
     * @return the gtgl parameter or {@code null} if no gtgl parameter was
     * specified.
     */
    public File gtgl() {
        return gtgl;
    }

    /**
     * Returns the ref parameter or {@code null} if no ref parameter was
     * specified.
     * @return the ref parameter or {@code null} if no ref parameter was
     * specified.
     */
    public File ref() {
        return ref;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter.
     */
    public String out() {
        return out;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     */
    public File excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    /**
     * Returns the ped parameter or {@code null}
     * if no ped parameter was specified.
     *
     * @return the ped parameter or {@code null}
     * if no ped parameter was specified.
     */
    public File ped() {
        return ped;
    }

    /**
     * Returns the chrom parameter or {@code null}
     * if no chrom parameter was specified.
     *
     * @return the chrom parameter or {@code null}
     * if no chrom parameter was specified.
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the maxlr parameter.
     * @return the maxlr parameter
     */
    public float maxlr() {
        return maxlr;
    }

    // algorithm parameters

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter.
     */
    public int nthreads() {
        return nthreads;
    }

    /**
     * Returns the window parameter.
     * @return the window parameter.
     */
    public int window() {
        return window;
    }

    /**
     * Return the overlap parameter.
     * @return the overlap parameter.
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the gprobs parameter.
     * @return the gprobs parameter.
     */
    public boolean gprobs() {
        return gprobs;
    }

    /**
     * Returns the impute parameter.
     * @return the impute parameter.
     */
    public boolean impute() {
        return impute;
    }

    /**
     * Returns the usephase parameter.
     * @return the usephase parameter.
     */
    public boolean usephase() {
        return usephase;
    }

    /**
     * Returns the singlescale parameter.
     * @return the singlescale parameter.
     */
    public float singlescale() {
        return singlescale;
    }

    /**
     * Returns the duocale parameter.
     * @return the duoscale parameter.
     */
    public float duoscale() {
        return duoscale;
    }

    /**
     * Returns the trioscale parameter.
     * @return the trioscale parameter.
     */
    public float trioscale() {
        return trioscale;
    }

    /**
     * Returns the burnin-its parameter.
     * @return the burnin-its parameter.
     */
    public int burnin_its() {
        return burnin_its;
    }

    /**
     * Returns the phase-its parameter.
     * @return the phase-its parameter.
     */
    public int phase_its() {
        return phase_its;
    }

    /**
     * Returns the impute-its parameter.
     * @return the impute-its parameter.
     */
    public int impute_its() {
        return impute_its;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter.
     */
    public long seed() {
        return seed;
    }

    // ibd parameters

    /**
     * Returns the ibd parameter.
     * @return the ibd parameter.
     */
    public boolean ibd() {
        return ibd;
    }

    /**
     * Returns the ibdlod parameter.
     * @return the ibdlod parameter
     */
    public float ibdlod() {
        return ibdlod;
    }

    /**
     * Returns the ibdscale parameter.
     * @return the ibdscale parameter.
     */
    public float ibdscale() {
        return ibdscale;
    }

    /**
     * Returns the ibdtrim parameter.
     * @return the ibdtrim parameter.
     */
    public int ibdtrim() {
        return ibdtrim;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter.
     */
    public File map() {
        return map;
    }

    // expert parameters

    /**
     * Return the nsamples parameter.
     * @return the nsamples parameter.
     */
    public int nsamples() {
        return nsamples;
    }

    /**
     * Returns the buildwindow parameter.
     * @return the buildwindow parameter.
     */
    public int buildwindow() {
        return buildwindow;
    }

    /**
     * Returns the simped parameter
     * @return the simped parameter
     */
    public boolean simped() { return simped; }

    /**
     * Returns the hapfrq parameter
     * @return the hapfrq parameter
     */
    public float hapfrq() { return hapfrq; }

    public File getFounders() { return founders; }

    public boolean jaccard() { return jaccard; }

    public boolean jacquard() { return jacquard; }

    public int advanceSize() { return advancesize; }

    public int skipSize() { return skipsize; }

    public int extendSize() { return extendsize; }

    public int minNHaps() { return minnhaps; }

    public float fittol() { return fittol; }

    public boolean constrain() { return constrain; }

    public boolean emitjmd() { return emitjmt; }
}
