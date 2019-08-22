package pedsim;

import beagleutil.BeagleUtils;
import beagleutil.Samples;
import dag.Dag;
import dag.MergeableDag;
import haplotype.*;
import main.BasicGenotypeValues;
import main.GenotypeValues;
import main.NuclearFamilies;
import main.Parameters;
import sample.TrioBaum;
import vcf.*;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

/**
 * Created by chartl on 7/2/2015.
 *
 * Simulates inherited alleles from a pedigree, with no constraints on recombination
 * (so allelic blocks are inherited independently, effectively allowing multiple
 * recombinations)
 *
 * The simulator takes the input *phased* reference VCF and a pedigree (for which the
 * reference samples are founders), and writes the founders and their (simulated) descendents
 * into the output VCF.
 *
 * The main parameters are the minimum marker frequency (hapfrq) and the block length. The minimum marker frequency
 * is useful because we don't necessarily want haplotypes to be tagged by single variants, but to simulate data
 * where they will have to be phased.
 *
 * Sampling from a pedigree is basically a run of the forward algorithm. The founders construct the DAG and have
 * perfect emission probabilities; children have uniformly uncertain genotypes at every marker.
 */
public class PedigreeSimulator {

    private Data data;
    private static PrintWriter debug;
    private static final boolean DEBUG = true;

    private PedigreeSimulator(Data data) {
        this.data = data;
    }

    protected static void makeSelection(short[][] trioSelect, Random random) {
        for ( int tId = 0; tId < trioSelect.length; tId ++ ) {
            trioSelect[tId][0] = (short) (random.nextBoolean() ? 0 : 1);
            trioSelect[tId][1] = (short) (random.nextBoolean() ? 0 : 1);
        }
    }

    protected static byte[] extractAlleles(HapPair sample, int which) {
        byte[] alleles = new byte[sample.nMarkers()];
        for ( int m = 0; m < alleles.length; m++ ) {
            if ( which == 0 ) {
                alleles[m] = sample.allele1(m);
            } else {
                alleles[m] = sample.allele2(m);
            }
        }
        return alleles;
    }

    private static void assert_(boolean bool) {
        if ( ! bool ) {
            throw new IllegalStateException("Assertion failure!");
        }
    }

    public static void runSimulation(Parameters params, Data data) {
        try {
            debug = new PrintWriter(new java.io.PrintStream("firstHaps.dbg"));
        } catch ( FileNotFoundException e ) {
            throw new IllegalStateException("Fix your shit");
        }
        AllData castData = (AllData) data;
        Samples allSamples = BeagleUtils.samplesFromPed(params.ped());
        NuclearFamilies nuclearFamilies = new NuclearFamilies(allSamples, params.ped());
        PrintWriter writer;
        try {
            writer = new PrintWriter(params.out());
        } catch (FileNotFoundException e) {
            writer = new PrintWriter(System.out);
        }
        VcfWriter.writeMetaLines(allSamples.ids(), "SimulatePedigree", true, true, false, writer);
        Random random = new Random(0xCAFECAFEL);  // to choose window sizes
        short[][] trioSelect = new short[nuclearFamilies.nTrios()][2];
        while (castData.canAdvanceRefWindow()) {
            castData.advanceRefWindow(0, 500);
            System.out.printf("at: %d:%d%n",
                    castData.refHaps().get(0).markers().markers()[0].chromIndex(),
                    castData.refHaps().get(0).markers().markers()[0].pos());
            HapPair[] genotypes = getGenotypes(allSamples, nuclearFamilies, castData, trioSelect);
            SampleHapPairs sampleHaplotypes = new BasicSampleHapPairs(allSamples, Arrays.asList(genotypes));
            GenotypeValues markerValues = new BasicGenotypeValues(sampleHaplotypes.markers(), allSamples);
            VcfWriter.appendRecords(sampleHaplotypes, markerValues, 0, sampleHaplotypes.nMarkers(), writer, false, allSamples.ids());
            makeSelection(trioSelect, random);
        }
        writer.close();

    }

    protected static HapPair[] getGenotypes(Samples all, NuclearFamilies families, AllData data, short[][] select) {
        HapPair[] genotypes = new HapPair[families.nSamples()];
        //System.out.printf("ref==%d, all==%d%n", data.refSamples().nSamples(), all.nSamples());
        int refIdx = 0;
        HashSet<String> refHaps = new HashSet<>();
        for ( HapPair hp : data.refHaps() ) {
            String refId = data.refSamples().id(refIdx++);
            int allIdx = all.index(refId);
            genotypes[allIdx] = hp;
            refHaps.addAll(Arrays.asList(h2s(hp).split("\\n")));
        }
        HashSet<String> hapCheck = new HashSet<String>();
        for ( int trioIdx = 0; trioIdx < families.nTrios(); trioIdx++ ) {
            HapPair dad = genotypes[all.index(families.trioFather(trioIdx))];
            HapPair mom = genotypes[all.index(families.trioMother(trioIdx))];
            byte[] dadAlleles = extractAlleles(dad, select[trioIdx][0]);
            byte[] momAlleles = extractAlleles(mom, select[trioIdx][1]);
            int childIdx = all.index(families.trioOffspring(trioIdx));
            HapPair child = new BitHapPair(dad.markers(), childIdx,
                    dadAlleles, momAlleles);
            genotypes[childIdx] = child;
            if ( DEBUG ) {
                String dadId = all.id(families.trioFather(dad.idIndex()));
                String momId = all.id(families.trioMother(mom.idIndex()));
                String childId = all.id(childIdx);
                debug.printf("Trio %d%n", trioIdx);
                debug.printf("Dad (%s):%n%s%n", dadId, h2s(dad));
                debug.printf("Mom (%s):%n%s%n", momId, h2s(mom));
                debug.printf("Child (%s):%n%s%n", childId, h2s(child)); // hmm
            }
            hapCheck.addAll(Arrays.asList(h2s(dad).split("\\n")));
            hapCheck.addAll(Arrays.asList(h2s(mom).split("\\n")));
            hapCheck.addAll(Arrays.asList(h2s(child).split("\\n")));
        }
        if ( hapCheck.size() != refHaps.size() ) {
            throw new IllegalStateException(String.format("Error: ref haps=%d, trio haps=%d", refHaps.size(), hapCheck.size()));
        }
        return genotypes;
    }

    private static String h2s(HapPair h) {
        StringBuilder sb1 = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        for ( int m = 0; m < h.nMarkers(); m ++ ) {
            sb1.append(h.allele1(m));
            sb2.append(h.allele2(m));
        }
        return String.format("%s%n%s", sb1.toString(), sb2.toString());
    }

}
