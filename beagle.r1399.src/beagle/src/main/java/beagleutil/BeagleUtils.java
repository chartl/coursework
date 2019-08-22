package beagleutil;

import blbutil.FileIterator;
import blbutil.InputIterator;
import blbutil.StringUtil;
import haplotype.HapPair;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Static class for utilities
 */
public class BeagleUtils {

    private BeagleUtils() {
        // no instantiation
    }


    public static Samples samplesFromPed(File pedFile) {
        List<String> sampleIds = new ArrayList<String>(100);
        try(FileIterator<String> pedIter = InputIterator.fromGzipFile(pedFile)) {
            while (pedIter.hasNext()) {
                String[] fields = getPedFields(pedIter.next().trim());
                sampleIds.add(fields[1]);
            }
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("Bad ped file");
        }
        String[] ids = new String[sampleIds.size()];
        sampleIds.toArray(ids);
        return Samples.fromIds(ids);
    }

    public static String[] getHaplotypes(HapPair pair) {
        StringBuilder b1 = new StringBuilder();
        StringBuilder b2 = new StringBuilder();
        for ( int m = 0; m < pair.nMarkers(); m ++ ) {
            b1.append(pair.allele1(m));
            b2.append(pair.allele2(m));
        }
        return new String[]{b1.toString(), b2.toString()};
    }

    public static int countUniqueHaplotypes(List<HapPair> hapPairs) {
        Set<String> haplos = new HashSet<>(2 * hapPairs.size());
        for (HapPair hap : hapPairs ) {
            StringBuilder b1 = new StringBuilder();
            StringBuilder b2 = new StringBuilder();
            for ( int m = 0; m < hap.nMarkers() ; m ++ ) {
                b1.append(hap.allele1(m));
                b2.append(hap.allele2(m));
            }
            haplos.add(b1.toString());
            haplos.add(b2.toString());
        }
        return haplos.size();
    }


    private static String[] getPedFields(String line) {
        String[] fields = StringUtil.getFields(line, 5);
        if (fields.length < 4) {
            String s = "invalid line in ped file: " + line;
            throw new IllegalArgumentException(s);
        }
        return fields;
    }
}
