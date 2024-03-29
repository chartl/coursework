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
package sample;

import dag.Dag;
import haplotype.HapPair;
import haplotype.BitHapPair;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import vcf.GL;

/**
 * Class {@code SingleBaum} implements the Baum forward and backward
 * algorithms for a hidden Markov model (HMM) of an individual's genotype data.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SingleBaum implements SingleBaumInterface {

    private final Dag dag;
    private final GL gl;
    private final int nMarkers;
    private final int nCopies;
    private final long seed;
    private final Random random;

    private final int[] node1;
    private final int[] node2;
    private final double[] nodeValue;

    private final byte[][] alleles1;
    private final byte[][] alleles2;

    private final SingleBaumLevel[] levels;
    private final SingleNodes fwdNodes;
    private final SingleNodes bwdNodes;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a new {@code SingleBaum} instance.
     *
     * @param dag the directed acyclic graph that determines the
     * transition probabilities.
     * @param gl the emission probabilities.
     * @param seed the initial random seed.
     * @param nCopies the number of haplotype pairs that will be sampled for
     * each individual.
     *
     * @throws IllegalArgumentException
     * if {@code nCopies<1 || dag.markers().equals(gl.markers())==false}
     * @throws NullPointerException if {@code dag==null || gl==null}
     */
    public SingleBaum(Dag dag, GL gl, long seed, int nCopies) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (nCopies < 1) {
            throw new IllegalArgumentException("nCopies<1: " + nCopies);
        }
        this.dag = dag;
        this.gl = gl;
        this.nMarkers = dag.nMarkers();
        this.nCopies = nCopies;
        this.seed = seed;
        this.random = new Random(seed);

        this.node1 = new int[nCopies];
        this.node2 = new int[nCopies];
        this.nodeValue = new double[nCopies];
        this.alleles1 = new byte[nCopies][gl.nMarkers()];
        this.alleles2 = new byte[nCopies][gl.nMarkers()];

        int size = (int) Math.ceil(Math.sqrt(1 + 8*dag.nMarkers())/2.0) + 1;
        this.levels = new SingleBaumLevel[size];
        for (int j=0; j<levels.length; ++j) {
            levels[j] = new SingleBaumLevel(dag, gl);
        }
        this.fwdNodes = new SingleNodes();
        this.bwdNodes = new SingleNodes();
    }

    @Override
    public Dag dag() {
        return dag;
    }

    @Override
    public GL gl() {
        return gl;
    }

    @Override
    public int nCopies() {
        return nCopies;
    }

    @Override
    public long seed() {
        return seed;
    }

    @Override
    public List<HapPair> randomSample(int sample) {
        forwardAlgorithm(sample);
        initSampleAlleles(currentLevel(), sample);
        for (int j=nMarkers-2; j>=0; --j) {
            SingleBaumLevel level = previousLevel(sample);
            sampleAlleles(level, sample);
        }
        return hapList(sample);
    }

    @Override
    public List<HapPair> randomSample(int sample, double[] gtProbs) {
        checkGtProbs(gtProbs);
        forwardAlgorithm(sample);
        initSampleAlleles(currentLevel(), sample);
        currentLevel().setInitialBackwardValues(bwdNodes);
        setGtProbs(currentLevel(), gtProbs);
        for (int j=nMarkers-2; j>=0; --j) {
            SingleBaumLevel level = previousLevel(sample);
            sampleAlleles(level, sample);
            level.setBackwardValues(bwdNodes);
            setGtProbs(level, gtProbs);
        }
        return hapList(sample);
    }

    private void checkGtProbs(double[] gtProbs) {
        if (gtProbs.length != gl.markers().sumGenotypes()) {
            String s = "gtProbs.length!=gl.markers().sumGenotypes()";
            throw new IllegalArgumentException(s);
        }
    }

    private void setGtProbs(SingleBaumLevel level, double[] gtProbs) {
        int m = level.marker();
        int nGenotypes = gl.marker(m).nGenotypes();
        int base = gl.markers().sumGenotypes(m);
        for (int j=0; j<nGenotypes; ++j) {
            gtProbs[base + j] = level.gtProbs(j);
        }
    }

    private List<HapPair> hapList(int sample) {
        List<HapPair> hapList = new ArrayList<>(2*nCopies);
        for (int copy=0; copy<nCopies; ++copy) {
            HapPair haps = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sample),
                    alleles1[copy], alleles2[copy]);
            hapList.add(haps);
        }
        return hapList;
    }

    private void initSampleAlleles(SingleBaumLevel level, int sample) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = initialRandomState(level);
            node1[copy] = level.parentNode1(state);
            node2[copy] = level.parentNode2(state);
            nodeValue[copy] =  parentSum(level, sample, state);
            alleles1[copy][m] = level.symbol1(state);
            alleles2[copy][m] = level.symbol2(state);
        }
    }

    private int initialRandomState(SingleBaumLevel level) {
        double d = random.nextDouble();
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            sum += level.forwardValue(j);
            if (d <= sum) {
                return j;
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private double parentSum(SingleBaumLevel level, int sample, int state) {
        int marker = level.marker();
        double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
        int edge1 = level.edge1(state);
        int edge2 = level.edge2(state);
        double tp1 = dag.condEdgeProb(marker, edge1);
        double tp2 = dag.condEdgeProb(marker, edge2);
        byte symbol1 = dag.symbol(marker, edge1);
        byte symbol2 = dag.symbol(marker, edge2);
        double ep = gl.gl(marker, sample, symbol1, symbol2);
        return fwdValue / ( ep*tp1*tp2 );
    }

    private void sampleAlleles(SingleBaumLevel level, int sample) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = randomPreviousState(level, node1[copy], node2[copy],
                    nodeValue[copy]);
            node1[copy] = level.parentNode1(state);
            node2[copy] = level.parentNode2(state);
            nodeValue[copy] =  parentSum(level, sample, state);
            alleles1[copy][m] = level.symbol1(state);
            alleles2[copy][m] = level.symbol2(state);
        }
    }

    private int randomPreviousState(SingleBaumLevel level, int node1,
            int node2, double nodeValue) {
        double d = random.nextDouble() * nodeValue;
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            if ( node1==level.childNode1(j)
                    && node2==level.childNode2(j) ) {
                sum += level.forwardValue(j);
                if (d <= sum) {
                    return j;
                }
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private SingleBaumLevel nextLevel() {
        ++arrayIndex;
        if (arrayIndex == levels.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return levels[arrayIndex];
    }

    private SingleBaumLevel currentLevel() {
        return levels[arrayIndex];
    }

    private SingleBaumLevel previousLevel(int sample) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            levels[arrayIndex].setChildNodes(fwdNodes);
            int startLevel = levels[windowIndex].marker() + 1;
            int endLevel = startLevel + (levels.length - (windowIndex + 1) );
            for (int marker=startLevel; marker<endLevel; ++marker) {
                nextLevel().setForwardValues(fwdNodes, marker, sample);
            }
            return currentLevel();
        }
        else {
            return levels[--arrayIndex];
        }
    }

    private void forwardAlgorithm(int sample) {
        SingleBaumLevel.initializeNodes(fwdNodes);
        this.windowIndex = -1;
        this.arrayIndex = levels.length - 1;
        for (int marker=0; marker<nMarkers; ++marker) {
            nextLevel().setForwardValues(fwdNodes, marker, sample);
        }
    }
}
