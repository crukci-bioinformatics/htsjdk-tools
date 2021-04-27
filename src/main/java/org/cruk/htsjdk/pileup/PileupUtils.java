/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk.pileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Utility methods for filtering records overlapping a locus and tabulating
 * pileup base counts.
 *
 * @author eldrid01
 */
public class PileupUtils {

    public static final byte[] VALID_BASES = new byte[] { 'A', 'C', 'G', 'T', 'N' };

    // lookup of bases to an index used in arrays of base counts
    private static final int BASES_ARRAY_LENGTH = 127;
    private static final int[] BASE_INDEX_LOOKUP = new int[BASES_ARRAY_LENGTH];
    static {
        Arrays.fill(BASE_INDEX_LOOKUP, -1);
        int shiftToLowerCase = 'a' - 'A';
        for (int i = 0; i < VALID_BASES.length; i++) {
            byte base = VALID_BASES[i];
            BASE_INDEX_LOOKUP[base] = i;
            BASE_INDEX_LOOKUP[base + shiftToLowerCase] = i;
        }
    }

    /**
     * Returns the index for accessing the count for the specific base in the array
     * returned by getBaseCounts.
     *
     * Valid bases are 'A', 'C', 'G', 'T', 'N' and their lower case equivalents.
     *
     * @param base the base
     * @return an index into the array of base counts
     */
    public static int getBaseIndex(byte base) {
        return BASE_INDEX_LOOKUP[base];
    }

    /**
     * Returns an array of base counts for the list of RecordAndOffset objects.
     *
     * @param pileup list of RecordAndOffset objects
     * @return an array of A, C, G, T and N counts
     */
    public static <T extends AbstractRecordAndOffset> int[] getBaseCounts(List<T> pileup) {
        int[] counts = new int[VALID_BASES.length];
        for (T recordAndOffset : pileup) {
            byte base = recordAndOffset.getReadBase();
            int index = BASE_INDEX_LOOKUP[base];
            if (index != -1) {
                counts[index]++;
            }
        }
        return counts;
    }

    /**
     * Filters a list of RecordAndOffset objects removing those with mapping quality
     * scores below the minimum value specified.
     *
     * @param pileup                list of RecordAndOffset objects
     * @param minimumMappingQuality minimum mapping quality score
     * @return filtered list of RecordAndOffset objects
     */
    public static <T extends AbstractRecordAndOffset> List<T> filterLowMappingQualities(List<T> pileup,
            int minimumMappingQuality) {
        List<T> filteredPileup = new ArrayList<>();
        for (T recordAndOffset : pileup) {
            SAMRecord record = recordAndOffset.getRecord();
            if (record.getMappingQuality() >= minimumMappingQuality) {
                filteredPileup.add(recordAndOffset);
            }
        }
        return filteredPileup;
    }

    /**
     * Filters a list of RecordAndOffset objects removing those with base quality
     * scores below the minimum value specified.
     *
     * @param pileup             list of RecordAndOffset objects
     * @param minimumBaseQuality minimum base quality score
     * @return filtered list of RecordAndOffset objects
     */
    public static <T extends AbstractRecordAndOffset> List<T> filterLowBaseQualities(List<T> pileup,
            int minimumBaseQuality) {
        List<T> filteredPileup = new ArrayList<>();
        for (T recordAndOffset : pileup) {
            if (recordAndOffset.getBaseQuality() >= minimumBaseQuality) {
                filteredPileup.add(recordAndOffset);
            }
        }
        return filteredPileup;
    }

    /**
     * Filters a list of RecordAndOffset objects removing those with base quality
     * scores and/or mapping qualities below the minimum values specified.
     *
     * @param pileup                list of RecordAndOffset objects
     * @param minimumBaseQuality    minimum base quality score
     * @param minimumMappingQuality minimum mapping quality score
     * @return filtered list of RecordAndOffset objects
     */
    public static <T extends AbstractRecordAndOffset> List<T> filterLowQualityScores(List<T> pileup,
            int minimumBaseQuality, int minimumMappingQuality) {
        List<T> filteredPileup = new ArrayList<>();
        for (T recordAndOffset : pileup) {
            SAMRecord record = recordAndOffset.getRecord();
            if (recordAndOffset.getBaseQuality() >= minimumBaseQuality
                    && record.getMappingQuality() >= minimumMappingQuality) {
                filteredPileup.add(recordAndOffset);
            }
        }
        return filteredPileup;
    }

    /**
     * Filters a list of RecordAndOffset objects for overlapping records for the
     * same read.
     *
     * At most, one record is retained for each read that overlaps this locus in
     * question that being the one with the highest base quality at that locus.
     * However, if records for the same read have different base calls at the locus,
     * all are discarded.
     *
     * @param pileup list of RecordAndOffset objects
     * @return filtered list of RecordAndOffset objects
     */
    public static <T extends AbstractRecordAndOffset> List<T> filterOverlaps(List<T> pileup) {

        Map<String, T> retained = new HashMap<>();
        Set<String> excluded = new HashSet<>();

        for (T recordAndOffset : pileup) {

            SAMRecord record = recordAndOffset.getRecord();
            String name = record.getReadName();

            // check if we've already seen records for this read that have
            // discordant base calls at this locus
            if (excluded.contains(name)) {
                continue;
            }

            // check if consistent with previously seen record for this read
            T previous = retained.get(name);
            if (previous == null) {
                retained.put(name, recordAndOffset);
            } else {
                if (SequenceUtil.basesEqual(recordAndOffset.getReadBase(), previous.getReadBase())) {
                    // retain record with highest base quality
                    if (recordAndOffset.getBaseQuality() > previous.getBaseQuality()) {
                        retained.put(name, recordAndOffset);
                    }
                } else {
                    retained.remove(name);
                    excluded.add(name);
                }
            }
        }

        return new ArrayList<>(retained.values());
    }

    /**
     * Filters a list of RecordAndOffset objects retaining only those for reads
     * belonging to read groups matching the given sample names.
     *
     * @param pileup  list of RecordAndOffset objects
     * @param samples the names of samples
     * @return filtered list of RecordAndOffset objects
     */
    public static <T extends AbstractRecordAndOffset> List<T> filterSamples(List<T> pileup, Set<String> samples) {
        List<T> filteredPileup = new ArrayList<>();
        for (T recordAndOffset : pileup) {
            SAMRecord record = recordAndOffset.getRecord();
            if (samples.contains(record.getReadGroup().getSample())) {
                filteredPileup.add(recordAndOffset);
            }
        }
        return filteredPileup;
    }
}
