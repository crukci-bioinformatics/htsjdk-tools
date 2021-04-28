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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Class representing the pileup of reads at a given locus as returned by e.g.
 * htsjdk.samtools.util.SamLocusIterator.LocusInfo.getRecordAndOffsets().
 *
 * No checks are made that the elements all refer to the same locus.
 *
 * @param <T> the base type of elements in the pileup
 */
public class Pileup<T extends AbstractRecordAndOffset> implements Iterable<T> {

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

    private List<T> elements = new ArrayList<>();

    /**
     * Creates a new empty Pileup object.
     */
    public Pileup() {
    }

    /**
     * Creates a new Pileup object populated with the given set of elements.
     *
     * @param elements
     */
    public Pileup(Collection<T> elements) {
        this.elements.addAll(elements);
    }

    /**
     * Adds an element to the pileup.
     *
     * @param element the pileup element.
     */
    public void addElement(T element) {
        elements.add(element);
    }

    /**
     * Returns the number of elements in the pileup.
     *
     * @return the pileup size
     */
    public int size() {
        return elements.size();
    }

    /**
     * Returns an array of base counts.
     *
     * @return an array of A, C, G, T and N counts
     */
    public int[] getBaseCounts() {
        int[] counts = new int[VALID_BASES.length];
        for (T element : elements) {
            byte base = element.getReadBase();
            int index = BASE_INDEX_LOOKUP[base];
            if (index != -1) {
                counts[index]++;
            }
        }
        return counts;
    }

    /**
     * Returns the number of elements having the specified base.
     *
     * @param base the base
     * @return the number of elements for the given base
     */
    public int getBaseCount(byte base) {
        int count = 0;
        for (T element : elements) {
            if (element.getReadBase() == base) {
                count++;
            }
        }
        return count;
    }

    /**
     * Returns a filtered pileup in which reads with a mapping quality below the
     * given threshold are removed.
     *
     * @param minimumMappingQuality the minimum mapping quality of retained reads
     * @return the filtered pileup
     */
    public Pileup<T> getMappingQualityFilteredPileup(int minimumMappingQuality) {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (element.getRecord().getMappingQuality() >= minimumMappingQuality) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    /**
     * Returns a filtered pileup in which reads with a base quality at the locus
     * below the given threshold are removed.
     *
     * @param minimumBaseQuality the minimum base quality of retained reads
     * @return the filtered pileup
     */
    public Pileup<T> getBaseQualityFilteredPileup(int minimumBaseQuality) {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (element.getBaseQuality() >= minimumBaseQuality) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    /**
     * Returns a filtered pileup in which reads with a mapping quality below the
     * given threshold or a base quality at the locus below the minimum value
     * specified are removed.
     *
     * @param minimumMappingQuality the minimum mapping quality of retained reads
     * @param minimumBaseQuality    the minimum base quality of retained reads
     * @return the filtered pileup
     */
    public Pileup<T> getMappingAndBaseQualityFilteredPileup(int minimumMappingQuality, int minimumBaseQuality) {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (element.getRecord().getMappingQuality() >= minimumMappingQuality
                    && element.getBaseQuality() >= minimumBaseQuality) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    /**
     * Returns a filtered pileup in which overlapping reads at the given locus are
     * removed.
     *
     * At most, one record is retained for each read that overlaps this locus in
     * question that being the one with the highest base quality at that locus.
     * However, if records for the same read have different base calls at the locus,
     * all/both are discarded.
     *
     * @return the filtered pileup
     */
    public Pileup<T> getOverlapFilteredPileup() {

        Map<String, T> retained = new HashMap<>();
        Set<String> excluded = new HashSet<>();

        for (T element : elements) {

            SAMRecord record = element.getRecord();
            String name = record.getReadName();

            // check if we've already seen records for this read that have
            // discordant base calls at this locus
            if (excluded.contains(name)) {
                continue;
            }

            // check if consistent with previously seen record for this read
            T previous = retained.get(name);
            if (previous == null) {
                retained.put(name, element);
            } else {
                if (SequenceUtil.basesEqual(element.getReadBase(), previous.getReadBase())) {
                    // retain record with highest base quality
                    if (element.getBaseQuality() > previous.getBaseQuality()) {
                        retained.put(name, element);
                    }
                } else {
                    retained.remove(name);
                    excluded.add(name);
                }
            }
        }

        return new Pileup<>(retained.values());
    }

    /**
     * Returns a filtered pileup in which only reads belonging to the given set of
     * samples are retained.
     *
     * @param samples the sample names
     * @return the filtered pileup
     */
    public Pileup<T> getSampleFilteredPileup(Set<String> samples) {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (samples.contains(element.getRecord().getReadGroup().getSample())) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    /**
     * Returns a filtered pileup containing reads that align to the negative (or
     * reverse) strand.
     *
     * @return the filtered pileup
     */
    public Pileup<T> getNegativeStrandPileup() {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (element.getRecord().getReadNegativeStrandFlag()) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    /**
     * Returns a filtered pileup containing reads that have the given base at the
     * locus.
     *
     * @param base the base
     * @return the filtered pileup
     */
    public Pileup<T> getBaseFilteredPileup(byte base) {
        Pileup<T> filteredPileup = new Pileup<T>();
        for (T element : elements) {
            if (element.getReadBase() == base) {
                filteredPileup.addElement(element);
            }
        }
        return filteredPileup;
    }

    @Override
    public Iterator<T> iterator() {
        return elements.iterator();
    }
}
