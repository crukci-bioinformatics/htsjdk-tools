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
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
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
     * Returns an array of base counts for the list of RecordAndOffset objects.
     *
     * @param pileup
     * @return an array of A, C, G, T and N counts
     */
    public static int[] getBaseCounts(List<RecordAndOffset> pileup) {
        int[] counts = new int[VALID_BASES.length];
        for (RecordAndOffset recordAndOffset : pileup) {
            byte base = recordAndOffset.getReadBase();
            int index = BASE_INDEX_LOOKUP[base];
            if (index != -1) {
                counts[index]++;
            }
        }
        return counts;
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
    public static List<RecordAndOffset> filterLowQualityScores(List<RecordAndOffset> pileup, int minimumBaseQuality,
            int minimumMappingQuality) {
        List<RecordAndOffset> filteredPileup = new ArrayList<RecordAndOffset>();
        for (RecordAndOffset recordAndOffset : pileup) {
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
    public static List<RecordAndOffset> filterOverlappingFragments(List<RecordAndOffset> pileup) {

        Map<String, RecordAndOffset> retained = new HashMap<String, RecordAndOffset>();
        Set<String> excluded = new HashSet<String>();

        for (RecordAndOffset recordAndOffset : pileup) {

            SAMRecord record = recordAndOffset.getRecord();
            String name = record.getReadName();

            // check if we've already seen records for this read that have
            // discordant base calls at this locus
            if (excluded.contains(name)) {
                continue;
            }

            // check if consistent with previously seen record for this read
            RecordAndOffset previous = retained.get(name);
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

        return new ArrayList<RecordAndOffset>(retained.values());
    }
}