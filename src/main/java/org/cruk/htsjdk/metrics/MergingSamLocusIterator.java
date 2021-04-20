/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk.metrics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.NoSuchElementException;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.LocusComparator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;

/**
 * Merging SAM locus iterator that combines LocusInfo data from multiple input
 * BAM files so that a single combined LocusInfo object is returned for each
 * locus.
 *
 * @author eldrid01
 */
public class MergingSamLocusIterator implements CloseableIterator<LocusInfo> {

    private SAMFileHeader header;
    private PeekableIterator<LocusInfo> iterator;

    /**
     * Creates a new merging iterator over the given loci that will produce
     * LocusInfo objects for each locus merged from multiple BAM files.
     *
     * @param bamFiles
     * @param intervals
     */
    public MergingSamLocusIterator(Collection<SamReader> readers, List<Interval> intervals) {

        // create merged header
        List<SAMFileHeader> headers = new ArrayList<>();
        for (SamReader reader : readers) {
            if (!reader.hasIndex()) {
                throw new SAMException("SAM readers must be indexed to create a merging SAM locus iterator");
            }
            headers.add(reader.getFileHeader());
        }
        SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SortOrder.coordinate, headers, true);
        header = headerMerger.getMergedHeader();

        // create sorted, unique interval list
        IntervalList intervalList = new IntervalList(header);
        intervalList.addall(intervals);
        intervalList = intervalList.uniqued();

        // create merging iterator
        List<CloseableIterator<LocusInfo>> locusIterators = new ArrayList<>();
        for (SamReader reader : readers) {
            SamLocusIterator locusIterator = new SamLocusIterator(reader, intervalList);
            locusIterators.add(locusIterator);
        }
        MergingIterator<LocusInfo> mergingLocusIterator = new MergingIterator<>(new LocusComparator<>(),
                locusIterators);
        iterator = new PeekableIterator<>(mergingLocusIterator);
    }

    @Override
    public boolean hasNext() {
        return iterator.hasNext();
    }

    @Override
    public LocusInfo next() {
        if (!iterator.hasNext())
            throw new NoSuchElementException();

        LocusInfo locusInfo = iterator.next();

        while (iterator.hasNext()) {
            LocusInfo nextLocusInfo = iterator.peek();
            if (!nextLocusInfo.getContig().equals(locusInfo.getContig())
                    || nextLocusInfo.getPosition() != locusInfo.getPosition()) {
                break;
            }
            addLocusInfo(locusInfo, nextLocusInfo);
            iterator.next();
        }

        return locusInfo;
    }

    @Override
    public void close() {
        iterator.close();
    }

    /**
     * @return the merged SAM header
     */
    public SAMFileHeader getHeader() {
        return header;
    }

    /**
     * Adds RecordAndOffset objects from one LocusInfo to another.
     *
     * @param locusInfo
     * @param additionalLocusInfo
     */
    private void addLocusInfo(LocusInfo locusInfo, LocusInfo additionalLocusInfo) {
        for (RecordAndOffset recordAndOffset : additionalLocusInfo.getRecordAndOffsets()) {
            locusInfo.add(recordAndOffset);
        }
        for (RecordAndOffset recordAndOffset : additionalLocusInfo.getDeletedInRecord()) {
            locusInfo.addDeleted(recordAndOffset.getRecord(), recordAndOffset.getOffset());
        }
        for (RecordAndOffset recordAndOffset : additionalLocusInfo.getInsertedInRecord()) {
            locusInfo.addInserted(recordAndOffset.getRecord(), recordAndOffset.getOffset());
        }
    }
}
