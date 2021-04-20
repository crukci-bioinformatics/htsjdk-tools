/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk.pileup;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.intervals.IntervalUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.util.SequenceUtil;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Command line tool for generating a pileup summary for all loci within a
 * specified set of intervals.
 *
 * @author eldrid01
 */
@Command(name = "pileup-counts", versionProvider = PileupCounts.class, description = "\nGenerates a pileup summary with read counts for each position and allele.\n", mixinStandardHelpOptions = true)
public class PileupCounts extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    @Option(names = "--id", description = "Identifier for this dataset; if included the pileup counts table will have an additional ID column (optional).")
    private String id;

    @Option(names = { "-i",
            "--input" }, required = true, description = "Input BAM file which must be in coordinate sort order and indexed (required).")
    private File bamFile;

    @Option(names = { "-l",
            "--intervals" }, required = true, description = "Intervals over which to generate pileup counts; can be in BED or Picard-style interval format (required).")
    private File intervalsFile;

    @Option(names = { "-r",
            "--reference-sequence" }, description = "Reference sequence FASTA file which must be indexed and have an accompanying dictionary (optional).")
    private File referenceSequenceFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "The output pileup summary table with read counts for each position and allele (required).")
    private File pileupCountsFile;

    @Option(names = { "-q",
            "--minimum-base-quality" }, description = "Minimum base quality for the base call at a locus for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumBaseQuality = 10;

    @Option(names = { "-m",
            "--minimum-mapping-quality" }, description = "Minimum mapping quality for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumMappingQuality = 1;

    @Option(names = "--output-N-counts", description = "Output counts for the number of N base calls.")
    private boolean outputNCounts = false;

    @Option(names = "--validation-stringency", description = "Validation stringency applied to the BAM file (default: ${DEFAULT-VALUE}).")
    private ValidationStringency validationStringency = ValidationStringency.LENIENT;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new PileupCounts()).execute(args);
        System.exit(exitCode);
    }

    /**
     * Main work method for iterating over loci within the given intervals and
     * tabulating read counts for each base.
     */
    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        ProgressLogger progress = new ProgressLogger(logger, 100, "loci");

        IOUtil.assertFileIsReadable(bamFile);
        IOUtil.assertFileIsReadable(intervalsFile);
        IOUtil.assertFileIsWritable(pileupCountsFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(validationStringency).open(bamFile);
        if (!reader.hasIndex()) {
            logger.error("No index found for input BAM file");
            return 1;
        }

        ReferenceSequenceFile referenceFile = null;
        if (referenceSequenceFile != null) {
            IOUtil.assertFileIsReadable(referenceSequenceFile);

            referenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequenceFile);

            if (!referenceFile.isIndexed()) {
                logger.error("Reference sequence FASTA file is not indexed");
                return 1;
            }

            // check sequence dictionary is consistent with the BAM file
            SAMSequenceDictionary dictionary = referenceFile.getSequenceDictionary();
            if (dictionary != null) {
                logger.info("Sequence dictionary found");
                if (!SequenceUtil.areSequenceDictionariesEqual(dictionary,
                        reader.getFileHeader().getSequenceDictionary())) {
                    if (validationStringency == ValidationStringency.LENIENT) {
                        logger.warn("Reference sequence dictionary not consistent with BAM file");
                    } else if (validationStringency == ValidationStringency.STRICT) {
                        logger.error("Reference sequence dictionary not consistent with BAM file");
                        return 1;
                    }
                }
            }
        }

        List<Interval> intervals = IntervalUtils.readIntervalFile(intervalsFile);
        IntervalList intervalList = new IntervalList(reader.getFileHeader());
        intervalList.addall(intervals);

        SamLocusIterator locusIterator = new SamLocusIterator(reader, intervalList);

        // exclude reads that are marked as failing platform/vendor quality checks
        locusIterator.setIncludeNonPfReads(false);

        // exclude secondary alignments and reads marked as duplicates
        locusIterator.setSamFilters(Arrays.asList(new SecondaryAlignmentFilter(), new DuplicateReadFilter()));

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(pileupCountsFile);

        writeHeader(writer);

        while (locusIterator.hasNext()) {
            LocusInfo locusInfo = locusIterator.next();

            String referenceBase = null;
            if (referenceFile != null) {
                referenceBase = referenceFile
                        .getSubsequenceAt(locusInfo.getContig(), locusInfo.getPosition(), locusInfo.getPosition())
                        .getBaseString();
            }

            List<RecordAndOffset> pileup = locusInfo.getRecordAndOffsets();

            List<RecordAndOffset> filteredPileup = PileupUtils.filterLowQualityScores(pileup, minimumBaseQuality,
                    minimumMappingQuality);

            filteredPileup = PileupUtils.filterOverlappingFragments(filteredPileup);

            writePileupCounts(writer, locusInfo, referenceBase, filteredPileup);

            progress.record(locusInfo.getContig(), locusInfo.getPosition());
        }

        writer.close();

        locusIterator.close();
        reader.close();
        if (referenceFile != null) {
            referenceFile.close();
        }

        logger.info("Finished");
        return 0;
    }

    /**
     * Write the header for the coverage/pileup output table.
     *
     * @param writer
     * @throws IOException
     */
    private void writeHeader(BufferedWriter writer) throws IOException {
        if (id != null) {
            writer.write("ID");
            writer.write("\t");
        }
        writer.write("Chromosome");
        writer.write("\t");
        writer.write("Position");
        if (referenceSequenceFile != null) {
            writer.write("\t");
            writer.write("Reference base");
        }
        writer.write("\t");
        writer.write("Depth unfiltered");
        writer.write("\t");
        writer.write("Depth");
        writer.write("\t");
        writer.write("A count");
        writer.write("\t");
        writer.write("C count");
        writer.write("\t");
        writer.write("G count");
        writer.write("\t");
        writer.write("T count");
        if (outputNCounts) {
            writer.write("\t");
            writer.write("N count");
        }
        writer.write("\n");
    }

    private void writePileupCounts(BufferedWriter writer, LocusInfo locusInfo, String referenceBase,
            List<RecordAndOffset> filteredPileup) throws IOException {

        if (id != null) {
            writer.write(id);
            writer.write("\t");
        }

        writer.write(locusInfo.getContig());
        writer.write("\t");
        writer.write(Integer.toString(locusInfo.getPosition()));

        if (referenceSequenceFile != null) {
            writer.write("\t");
            writer.write(referenceBase);
        }

        writer.write("\t");
        writer.write(Integer.toString(locusInfo.size()));

        writer.write("\t");
        writer.write(Integer.toString(filteredPileup.size()));

        int[] baseCounts = PileupUtils.getBaseCounts(filteredPileup);
        int n = baseCounts.length;
        if (!outputNCounts) {
            n -= 1;
        }
        for (int i = 0; i < n; i++) {
            writer.write("\t");
            writer.write(Integer.toString(baseCounts[i]));
        }

        writer.write("\n");
    }
}
