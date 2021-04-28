/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk.metrics;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.pileup.Pileup;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Command line tool for computing a set of metrics for single nucleotide
 * variants (SNVs) from a given VCF file. These metrics are added to the INFO
 * column in the output VCF file and can be used for subsequent filtering to
 * obtain a higher precision set of variants.
 *
 * @author eldrid01
 */
@Command(name = "calculate-snv-metrics", versionProvider = CalculateSnvMetrics.class, description = "\nCalculates metrics for single nucleotide variants (SNVs) from a VCF file based on aligned reads from the given BAM file(s).\n", mixinStandardHelpOptions = true)
public class CalculateSnvMetrics extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger();

    // metric names
    private static final String DEPTH = "Depth";
    private static final String DEPTH_CONTROL = "DepthControl";
    private static final String LOW_MAPPING_QUALITY = "LowMapQual";
    private static final String READ_COUNT = "ReadCount";
    private static final String READ_COUNT_CONTROL = "ReadCountControl";
    private static final String VARIANT_ALLELE_COUNT = "VariantAlleleCount";
    private static final String VARIANT_ALLELE_COUNT_CONTROL = "VariantAlleleCountControl";
    private static final String VARIANT_ALLELE_FREQUENCY = "VariantAlleleFrequency";
    private static final String VARIANT_ALLELE_FREQUENCY_CONTROL = "VariantAlleleFrequencyControl";
    private static final String STRAND_BIAS = "StrandBias";
    private static final String VARIANT_STRAND_BIAS = "VariantStrandBias";
    private static final String REFERENCE_STRAND_BIAS = "ReferenceStrandBias";
    private static final String VARIANT_BASE_QUALITY_MEAN = "VariantBaseQual";
    private static final String VARIANT_BASE_QUALITY_MEDIAN = "VariantBaseQualMedian";
    private static final String VARIANT_MAPPING_QUALITY_MEAN = "VariantMapQual";
    private static final String VARIANT_MAPPING_QUALITY_MEDIAN = "VariantMapQualMedian";
    private static final String MAPPING_QUALITY_MEAN_DIFFERENCE = "MapQualDiff";
    private static final String MAPPING_QUALITY_MEDIAN_DIFFERENCE = "MapQualDiffMedian";
    private static final String VARIANT_MMQS_MEAN = "VariantMMQS";
    private static final String VARIANT_MMQS_MEDIAN = "VariantMMQSMedian";
    private static final String MMQS_MEAN_DIFFERENCE = "MMQSDiff";
    private static final String MMQS_MEDIAN_DIFFERENCE = "MMQSDiffMedian";
    private static final String DISTANCE_TO_ALIGNMENT_END_MEAN = "DistanceToAlignmentEnd";
    private static final String DISTANCE_TO_ALIGNMENT_END_MEDIAN = "DistanceToAlignmentEndMedian";
    private static final String DISTANCE_TO_ALIGNMENT_END_MAD = "DistanceToAlignmentEndMAD";
    private static final String HOMOPOLYMER_LENGTH = "HomopolymerLength";
    private static final String REPEAT_LENGTH = "Repeat";

    // upper limits for repeat and homopolymer lengths to detect
    private static final int MAX_HOMOPOLYMER_LENGTH = 50;
    private static final int MAX_REPEAT_LENGTH = 30;

    @Option(names = { "-i",
            "--input" }, required = true, description = "Input BAM file(s) which must be in coordinate sort order and indexed (required).")
    private List<File> bamFiles;

    @Option(names = { "-v",
            "--variants" }, required = true, description = "Input VCF file containing single nucleotide variants (required).")
    private File inputVcfFile;

    @Option(names = { "-r",
            "--reference-sequence" }, required = true, description = "Reference sequence FASTA file which must be indexed and have an accompanying dictionary (required).")
    private File referenceSequenceFile;

    @Option(names = { "-o",
            "--output" }, required = true, description = "Output VCF file for which SNV metrics have been added as INFO fields (required).")
    private File outputVcfFile;

    @Option(names = "--sample", description = ".The sample(s) for which to compute SNV metrics. Can be specified multiple times for multiple samples; used to restrict the reads used in computing metrics, if not specified metrics are computed from all reads within the given BAM files (optional).")
    private HashSet<String> samples;

    @Option(names = "--control_sample", description = "The sample(s) for which to compute SNV metrics for the reference allele. Can be specified multiple times for multiple samples; if not specified metrics will be computed for all reference-supporting reads in the samples specified using the sample argument (optional).")
    private HashSet<String> controlSamples;

    @Option(names = { "-q",
            "--minimum-base-quality" }, description = "Minimum base quality for the base call at a locus for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumBaseQuality = 10;

    @Option(names = { "-m",
            "--minimum-mapping-quality" }, description = "Minimum mapping quality for reads to be included (default: ${DEFAULT-VALUE}).")
    private int minimumMappingQuality = 1;

    @Option(names = "--validation-stringency", description = "Validation stringency applied to the BAM file(s) (default: ${DEFAULT-VALUE}).")
    private ValidationStringency validationStringency = ValidationStringency.LENIENT;

    private ReferenceSequenceFile referenceSequence;

    public static void main(String[] args) {
        int exitCode = new CommandLine(new CalculateSnvMetrics()).execute(args);
        System.exit(exitCode);
    }

    /**
     * Main run method for iterating over loci within the given intervals and
     * tabulating read counts for each base.
     */
    @Override
    public Integer call() throws Exception {
        logger.info(getClass().getName() + " (" + getPackageNameAndVersion() + ")");

        ProgressLogger progress = new ProgressLogger(logger, 100, "loci");

        checkInputsAndOutputs();

        // create SAM readers
        List<SamReader> samReaders = new ArrayList<>();
        for (File bamFile : bamFiles) {
            SamReader reader = SamReaderFactory.makeDefault().validationStringency(validationStringency).open(bamFile);
            if (!reader.hasIndex()) {
                logger.error("No index found for input BAM file " + bamFile.getName());
                return 1;
            }
            samReaders.add(reader);
        }

        // read single nucleotide variants
        List<VariantContext> snvs = readSNVs();

        // create intervals
        List<Interval> snvLoci = createIntervals(snvs);

        // exclude supplementary and secondary alignments and reads marked as duplicates
        List<SamRecordFilter> samFilters = Arrays.asList(new SecondaryOrSupplementaryFilter(),
                new DuplicateReadFilter());

        // create merged SAM locus iterator
        MergingSamLocusIterator mergingSamLocusIterator = new MergingSamLocusIterator(samReaders, samFilters, snvLoci);

        SAMSequenceDictionary sequenceDictionary = mergingSamLocusIterator.getHeader().getSequenceDictionary();

        // create SNV iterator
        VariantContextComparator snvComparator = new VariantContextComparator(sequenceDictionary);
        snvs.sort(snvComparator);
        PeekableIterator<VariantContext> snvIterator = new PeekableIterator<>(snvs.iterator());

        // reference sequence
        referenceSequence = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequenceFile);
        if (!referenceSequence.isIndexed()) {
            logger.error("Reference sequence FASTA file is not indexed");
            return 1;
        }

        // check sequence dictionary is consistent with the BAM file(s)
        SAMSequenceDictionary dictionary = referenceSequence.getSequenceDictionary();
        if (dictionary != null) {
            if (!SequenceUtil.areSequenceDictionariesEqual(dictionary, sequenceDictionary)) {
                if (validationStringency == ValidationStringency.LENIENT) {
                    logger.warn("Reference sequence dictionary not consistent with BAM file");
                } else if (validationStringency == ValidationStringency.STRICT) {
                    logger.error("Reference sequence dictionary not consistent with BAM file");
                    return 1;
                }
            }
        }

        // metrics lookup
        Map<String, Map<String, Object>> metricsLookup = new HashMap<>();

        while (mergingSamLocusIterator.hasNext()) {
            LocusInfo locusInfo = mergingSamLocusIterator.next();

            String contig = locusInfo.getContig();
            int position = locusInfo.getPosition();
            byte[] referenceBases = referenceSequence.getSubsequenceAt(contig, position, position).getBases();

            Pileup<RecordAndOffset> pileup = new Pileup<>(locusInfo.getRecordAndOffsets());

            boolean foundMatchingSnv = false;

            while (isVariantAtLocus(snvIterator.peek(), locusInfo)) {
                foundMatchingSnv = true;

                VariantContext snv = snvIterator.next();

                if (!snv.getReference().basesMatch(referenceBases)) {
                    logger.error("Mismatching reference base for " + getSnvId(snv));
                    return 1;
                }

                String id = getSnvId(snv);

                // shouldn't really have multiple entries for the same variant
                // but only compute the metrics once for each distinct SNV
                if (!metricsLookup.containsKey(id)) {
                    Map<String, Object> metrics = calculateSnvMetrics(snv, pileup);
                    metricsLookup.put(id, metrics);
                }
            }

            if (!foundMatchingSnv) {
                logger.error("Unexpected mismatch between locus and SNV iterators");
                return 1;
            }

            progress.record(locusInfo.getContig(), locusInfo.getPosition());
        }

        snvIterator.close();
        mergingSamLocusIterator.close();
        for (SamReader reader : samReaders) {
            reader.close();
        }
        if (referenceSequence != null) {
            referenceSequence.close();
        }

        logger.info("Writing VCF file containing SNV metrics");
        writeVcf(metricsLookup);

        logger.info("Finished");
        return 0;
    }

    /**
     * Check input files are readable and outputs are writable.
     */
    private void checkInputsAndOutputs() {
        for (File bamFile : bamFiles) {
            IOUtil.assertFileIsReadable(bamFile);
        }
        IOUtil.assertFileIsReadable(inputVcfFile);
        IOUtil.assertFileIsReadable(referenceSequenceFile);
        IOUtil.assertFileIsWritable(outputVcfFile);
    }

    /**
     * Read single nucleotide variants (SNVs) from the input VCF file.
     *
     * @return a list of single nucleotide variants (SNVs)
     */
    private List<VariantContext> readSNVs() {
        List<VariantContext> snvs = new ArrayList<>();
        VCFFileReader reader = new VCFFileReader(inputVcfFile, false);
        for (VariantContext variant : reader) {
            if (variant.isSNP()) {
                snvs.add(variant);
                if (variant.getAlternateAlleles().size() > 1) {
                    logger.warn("Warning multiple alternate alleles for variant " + getSnvId(variant)
                            + "; only computing metrics for allele with highest allele count.");
                }
            }
        }
        reader.close();
        return snvs;
    }

    /**
     * Returns a list of intervals for each variant.
     *
     * @param variants
     * @return
     */
    private List<Interval> createIntervals(List<VariantContext> variants) {
        List<Interval> intervals = new ArrayList<>();
        for (VariantContext variant : variants) {
            intervals.add(new Interval(variant));
        }
        return intervals;
    }

    /**
     * @param variant
     * @param locus
     * @return true if the variant is at the given locus
     */
    private boolean isVariantAtLocus(VariantContext variant, LocusInfo locus) {
        return variant != null && variant.getContig().equals(locus.getContig())
                && variant.getStart() == locus.getPosition();
    }

    /**
     * Returns an identifier for the SNV that can be used as a key in a lookup.
     *
     * @param snv
     * @return
     */
    private String getSnvId(VariantContext snv) {
        return String.format("%s:%d %c>%c", snv.getContig(), snv.getStart(), snv.getReference().getBases()[0],
                snv.getAltAlleleWithHighestAlleleCount().getBases()[0]);
    }

    /**
     * Calculates SNV metrics for the given variant and read pileup.
     *
     * For multi-allelic variants, metrics are only calculated for the alternate
     * allele with the highest allele count.
     *
     * @param snv                   the SNV variant
     * @param pileup                the read pileup
     * @param referenceSequenceFile the reference sequence file
     * @return the computed SNV metrics
     */
    private Map<String, Object> calculateSnvMetrics(VariantContext snv, Pileup<RecordAndOffset> pileup) {

        String chromosome = snv.getContig();
        int position = snv.getStart();

        Map<String, Object> metrics = new HashMap<>();

        byte referenceBase = snv.getReference().getBases()[0];
        byte variantAllele = snv.getAltAlleleWithHighestAlleleCount().getBases()[0];

        // depth metrics calculated before applying additional mapping quality and
        // overlap filters
        metrics.put(DEPTH, samples.isEmpty() ? pileup.size() : pileup.getSampleFilteredPileup(samples).size());

        if (!controlSamples.isEmpty()) {
            metrics.put(DEPTH_CONTROL, pileup.getSampleFilteredPileup(controlSamples).size());
        }

        int unfilteredReadCount = pileup.size();

        // apply mapping quality filter
        pileup = pileup.getMappingQualityFilteredPileup(minimumMappingQuality);

        if (unfilteredReadCount > 0) {
            metrics.put(LOW_MAPPING_QUALITY, 1.0 - (pileup.size() / (double) unfilteredReadCount));
        }

        // apply base quality filter
        pileup = pileup.getBaseQualityFilteredPileup(minimumBaseQuality);

        // apply filter for overlapping fragments
        pileup = pileup.getOverlapFilteredPileup();

        // apply filter for the specified sample(s)
        Pileup<RecordAndOffset> samplePileup = samples.isEmpty() ? pileup : pileup.getSampleFilteredPileup(samples);

        // collect reference metrics from control samples if these exist,
        // otherwise from reference-supporting reads in the same set of
        // samples used for computing the variant metrics
        Pileup<RecordAndOffset> referencePileup = controlSamples.isEmpty() ? samplePileup
                : pileup.getSampleFilteredPileup(controlSamples);

        int readCount = samplePileup.size();
        metrics.put(READ_COUNT, readCount);

        int referenceReadCount = referencePileup.size();
        if (!controlSamples.isEmpty()) {
            metrics.put(READ_COUNT_CONTROL, referenceReadCount);
        }

        int variantReadCount = samplePileup.getBaseCount(variantAllele);
        metrics.put(VARIANT_ALLELE_COUNT, variantReadCount);
        if (readCount > 0) {
            metrics.put(VARIANT_ALLELE_FREQUENCY, variantReadCount / (double) readCount);
        }

        if (!controlSamples.isEmpty()) {
            int variantAlleleCountControl = referencePileup.getBaseCount(variantAllele);
            metrics.put(VARIANT_ALLELE_COUNT_CONTROL, variantAlleleCountControl);
            if (referenceReadCount > 0) {
                metrics.put(VARIANT_ALLELE_FREQUENCY_CONTROL, variantAlleleCountControl / (double) referenceReadCount);
            }
        }

        calculateStrandBiasMetrics(variantAllele, referenceBase, samplePileup, referencePileup, metrics);

        Pileup<RecordAndOffset> variantReadPileup = samplePileup.getBaseFilteredPileup(variantAllele);
        Pileup<RecordAndOffset> referenceReadPileup = referencePileup.getBaseFilteredPileup(referenceBase);

        calculateBaseQualityMetrics(variantReadPileup, metrics);
        calculateMismatchBaseQualityMetrics(variantReadPileup, referenceReadPileup, metrics);
        calculateMappingQualityMetrics(variantReadPileup, referenceReadPileup, metrics);
        calculateDistanceToAlignmentEndMetrics(variantReadPileup, metrics);

        metrics.put(HOMOPOLYMER_LENGTH, getHomopolymerLength(chromosome, position));
        metrics.put(REPEAT_LENGTH, getRepeatLength(chromosome, position));

        return metrics;
    }

    /**
     * Calculate strand bias metrics.
     *
     * @param variantAllele
     * @param referenceBase
     * @param samplePileup
     * @param referencePileup
     * @param metrics
     */
    private void calculateStrandBiasMetrics(byte variantAllele, byte referenceBase,
            Pileup<RecordAndOffset> samplePileup, Pileup<RecordAndOffset> referencePileup,
            Map<String, Object> metrics) {

        Pileup<RecordAndOffset> negativeStrandPileup = samplePileup.getNegativeStrandPileup();

        int sampleDepth = samplePileup.size();
        if (sampleDepth > 0) {
            double strandBias = negativeStrandPileup.size() / (double) sampleDepth;
            if (strandBias > 0.5) {
                strandBias = 1.0 - strandBias;
            }
            metrics.put(STRAND_BIAS, strandBias);
        }

        int variantReadCount = samplePileup.getBaseCount(variantAllele);
        if (variantReadCount > 0) {
            double variantStrandBias = negativeStrandPileup.getBaseCount(variantAllele) / (double) variantReadCount;
            if (variantStrandBias > 0.5) {
                variantStrandBias = 1.0 - variantStrandBias;
            }
            metrics.put(VARIANT_STRAND_BIAS, variantStrandBias);
        }

        int referenceDepth = referencePileup.size();
        if (referenceDepth > 0) {
            Pileup<RecordAndOffset> referenceNegativeStrandPileup = referencePileup.getNegativeStrandPileup();
            double referenceStrandBias = referenceNegativeStrandPileup.getBaseCount(referenceBase)
                    / (double) referenceDepth;
            if (referenceStrandBias > 0.5) {
                referenceStrandBias = 1.0 - referenceStrandBias;
            }
            metrics.put(REFERENCE_STRAND_BIAS, referenceStrandBias);
        }
    }

    /**
     * Calculates metrics based on base quality scores in the variant reads.
     *
     * @param variantReadPileup
     * @param metrics
     */
    private void calculateBaseQualityMetrics(Pileup<RecordAndOffset> variantReadPileup, Map<String, Object> metrics) {

        DescriptiveStatistics variantBaseQualityStats = new DescriptiveStatistics();

        for (RecordAndOffset recordAndOffset : variantReadPileup) {
            variantBaseQualityStats.addValue(recordAndOffset.getBaseQuality());
        }

        if (variantBaseQualityStats.getN() > 0) {
            metrics.put(VARIANT_BASE_QUALITY_MEAN, variantBaseQualityStats.getMean());
            metrics.put(VARIANT_BASE_QUALITY_MEDIAN, variantBaseQualityStats.getPercentile(50));
        }
    }

    /**
     * Calculates mismatch base quality score metrics.
     * 
     * @param variantReadPileup
     * @param referenceReadPileup
     * @param metrics
     */
    private void calculateMismatchBaseQualityMetrics(Pileup<RecordAndOffset> variantReadPileup,
            Pileup<RecordAndOffset> referenceReadPileup, Map<String, Object> metrics) {

        DescriptiveStatistics variantMMQSStatistics = new DescriptiveStatistics();
        DescriptiveStatistics referenceMMQSStats = new DescriptiveStatistics();

        for (RecordAndOffset recordAndOffset : variantReadPileup) {
            variantMMQSStatistics.addValue(getMismatchQualitySum(recordAndOffset));
        }

        for (RecordAndOffset recordAndOffset : referenceReadPileup) {
            referenceMMQSStats.addValue(getMismatchQualitySum(recordAndOffset));
        }

        if (variantMMQSStatistics.getN() > 0) {
            metrics.put(VARIANT_MMQS_MEAN, variantMMQSStatistics.getMean());
            metrics.put(VARIANT_MMQS_MEDIAN, variantMMQSStatistics.getPercentile(50));
            if (referenceMMQSStats.getN() > 0) {
                metrics.put(MMQS_MEAN_DIFFERENCE, variantMMQSStatistics.getMean() - referenceMMQSStats.getMean());
                metrics.put(MMQS_MEDIAN_DIFFERENCE,
                        variantMMQSStatistics.getPercentile(50) - referenceMMQSStats.getPercentile(50));
            }
        }
    }

    /**
     * Calculates metrics based on mapping qualities in the variant and reference
     * reads.
     *
     * @param variantReadPileup
     * @param referenceReadPileup
     * @param metrics
     */
    private void calculateMappingQualityMetrics(Pileup<RecordAndOffset> variantReadPileup,
            Pileup<RecordAndOffset> referenceReadPileup, Map<String, Object> metrics) {

        DescriptiveStatistics variantMappingQualityStats = new DescriptiveStatistics();
        DescriptiveStatistics referenceMappingQualityStats = new DescriptiveStatistics();

        for (RecordAndOffset recordAndOffset : variantReadPileup) {
            variantMappingQualityStats.addValue(recordAndOffset.getRecord().getMappingQuality());
        }

        for (RecordAndOffset recordAndOffset : referenceReadPileup) {
            referenceMappingQualityStats.addValue(recordAndOffset.getRecord().getMappingQuality());
        }

        if (variantMappingQualityStats.getN() > 0) {
            metrics.put(VARIANT_MAPPING_QUALITY_MEAN, variantMappingQualityStats.getMean());
            metrics.put(VARIANT_MAPPING_QUALITY_MEDIAN, variantMappingQualityStats.getPercentile(50));
            if (referenceMappingQualityStats.getN() > 0) {
                metrics.put(MAPPING_QUALITY_MEAN_DIFFERENCE,
                        referenceMappingQualityStats.getMean() - variantMappingQualityStats.getMean());
                metrics.put(MAPPING_QUALITY_MEDIAN_DIFFERENCE,
                        referenceMappingQualityStats.getPercentile(50) - variantMappingQualityStats.getPercentile(50));
            }
        }
    }

    /**
     * Calculates distance to alignment end metrics for the variant position within
     * reads.
     *
     * @param variantReadPileup
     * @param metrics
     */
    private void calculateDistanceToAlignmentEndMetrics(Pileup<RecordAndOffset> variantReadPileup,
            Map<String, Object> metrics) {

        DescriptiveStatistics variantDistanceToAlignmentEndStats = new DescriptiveStatistics();

        for (RecordAndOffset recordAndOffset : variantReadPileup) {
            variantDistanceToAlignmentEndStats.addValue(getDistanceToAlignmentEnd(recordAndOffset));
        }

        if (variantDistanceToAlignmentEndStats.getN() > 0) {
            metrics.put(DISTANCE_TO_ALIGNMENT_END_MEAN, variantDistanceToAlignmentEndStats.getMean());
            metrics.put(DISTANCE_TO_ALIGNMENT_END_MEDIAN, variantDistanceToAlignmentEndStats.getPercentile(50));
            metrics.put(DISTANCE_TO_ALIGNMENT_END_MAD, getMedianAbsoluteDeviation(variantDistanceToAlignmentEndStats));
        }
    }

    /**
     * Gets the mismatch quality sum for the given read alignment and offset.
     *
     * The mismatch quality sum is the sum of base qualities for positions in the
     * read that align with a mismatch to the reference sequence excluding the
     * offset position at which there is a possible variant.
     *
     * @param recordAndOffset
     * @return
     */
    private int getMismatchQualitySum(RecordAndOffset recordAndOffset) {

        int offset = recordAndOffset.getOffset();

        SAMRecord record = recordAndOffset.getRecord();

        byte[] readBases = record.getReadBases();
        byte[] baseQualities = record.getBaseQualities();

        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();

        byte[] referenceBases = referenceSequence
                .getSubsequenceAt(record.getReferenceName(), alignmentStart, alignmentEnd).getBases();

        int mmqs = 0;

        for (AlignmentBlock block : record.getAlignmentBlocks()) {
            int readBlockStart = block.getReadStart() - 1;
            int referenceStart = block.getReferenceStart();
            int referenceBlockStart = referenceStart - alignmentStart;
            for (int i = 0; i < block.getLength(); i++) {
                if ((readBlockStart + i) == offset) {
                    continue;
                }
                byte referenceBase = referenceBases[referenceBlockStart + i];
                byte readBase = readBases[readBlockStart + i];
                if (!SequenceUtil.isNoCall(readBase) && !SequenceUtil.isNoCall(referenceBase)
                        && !SequenceUtil.basesEqual(readBase, referenceBase)) {
                    mmqs += baseQualities[readBlockStart + i];
                }
            }
        }

        return mmqs;
    }

    /**
     * Computes the distance between the variant position with a read and the
     * closest end of the aligned portion of the read.
     *
     * The distance is 1-based in that the smallest distance will be 1 where the
     * position for the pileup element is either the first or the last aligned
     * position.
     *
     * @param recordAndOffset
     * @return
     */
    private int getDistanceToAlignmentEnd(RecordAndOffset recordAndOffset) {
        int positionInRead = recordAndOffset.getOffset() + 1;

        SAMRecord record = recordAndOffset.getRecord();

        List<AlignmentBlock> alignmentBlocks = record.getAlignmentBlocks();

        AlignmentBlock firstAlignmentBlock = alignmentBlocks.get(0);
        int alignedStartPositionInRead = firstAlignmentBlock.getReadStart();

        AlignmentBlock lastAlignmentBlock = alignmentBlocks.get(alignmentBlocks.size() - 1);
        int alignedEndPositionInRead = lastAlignmentBlock.getReadStart() + lastAlignmentBlock.getLength() - 1;

        return Math.min(positionInRead - alignedStartPositionInRead, alignedEndPositionInRead - positionInRead) + 1;
    }

    /**
     * Returns the length of the longest homopolymer adjacent to or surrounding the
     * given position.
     *
     * @param chromosome
     * @param position
     * @return
     */
    private int getHomopolymerLength(String chromosome, long position) {
        int referenceLength = referenceSequence.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        byte[] referenceBases = referenceSequence
                .getSubsequenceAt(chromosome, position, Math.min(referenceLength, position + MAX_HOMOPOLYMER_LENGTH))
                .getBases();
        byte referenceBase = referenceBases[0];
        byte rightBase = 0;
        int rightCount = 0;
        if (referenceBases.length > 1) {
            rightBase = referenceBases[1];
            rightCount = 1;
            for (int i = 2; i < referenceBases.length; i++) {
                if (!SequenceUtil.basesEqual(rightBase, referenceBases[i]))
                    break;
                rightCount++;
            }
            if (SequenceUtil.basesEqual(rightBase, referenceBase))
                rightCount++;
        }

        referenceBases = referenceSequence
                .getSubsequenceAt(chromosome, Math.max(1, position - MAX_HOMOPOLYMER_LENGTH), position).getBases();
        byte leftBase = 0;
        int leftCount = 0;
        if (referenceBases.length > 1) {
            leftBase = referenceBases[referenceBases.length - 2];
            leftCount = 1;
            for (int i = referenceBases.length - 3; i >= 0; i--) {
                if (!SequenceUtil.basesEqual(leftBase, referenceBases[i]))
                    break;
                leftCount++;
            }
            if (SequenceUtil.basesEqual(leftBase, referenceBase))
                leftCount++;
        }

        if (leftCount > 0 && rightCount > 0 && SequenceUtil.basesEqual(leftBase, referenceBase)
                && SequenceUtil.basesEqual(rightBase, referenceBase)) {
            return leftCount + rightCount - 1;
        } else {
            return Math.max(leftCount, rightCount);
        }
    }

    /**
     * Returns the length of repetitive sequence adjacent to the variant position
     * where repeats can be 1-, 2-, 3- or 4-mers.
     *
     * @param chromosome
     * @param position
     * @return
     */
    private int getRepeatLength(String chromosome, long position) {
        int repeat1 = getRepeatCount(1, chromosome, position);
        int repeat2 = getRepeatCount(2, chromosome, position);
        int repeat3 = getRepeatCount(3, chromosome, position);
        int repeat4 = getRepeatCount(4, chromosome, position);

        int repeat = 0;
        if (repeat1 > 1 || repeat2 > 1 || repeat3 > 1 || repeat4 > 1) {
            repeat = Math.max(repeat1, repeat2 * 2);
            repeat = Math.max(repeat, repeat3 * 3);
            repeat = Math.max(repeat, repeat4 * 4);
        }

        return repeat;
    }

    /**
     * Returns the maximum number of repeats of the specified length adjacent to the
     * given position.
     *
     * @param n          the repeat length, e.g. n = 2 for dinucleotide repeats
     * @param chromosome
     * @param position
     * @return
     */
    private int getRepeatCount(int n, String chromosome, long position) {
        int referenceLength = referenceSequence.getSequenceDictionary().getSequence(chromosome).getSequenceLength();

        int count = 0;

        if (position < referenceLength) {
            byte[] referenceBases = referenceSequence
                    .getSubsequenceAt(chromosome, position + 1, Math.min(referenceLength, position + MAX_REPEAT_LENGTH))
                    .getBases();
            count = Math.max(count, getRepeatCount(n, referenceBases));
        }

        if (position > 1) {
            byte[] referenceBases = referenceSequence
                    .getSubsequenceAt(chromosome, Math.max(1, position - MAX_REPEAT_LENGTH), position - 1).getBases();
            ArrayUtils.reverse(referenceBases);
            count = Math.max(count, getRepeatCount(n, referenceBases));
        }

        return count;
    }

    /**
     * Returns the number of repeats of length n in the given sequence.
     *
     * @param n     the repeat length, e.g. n = 2 for dinucleotide repeats
     * @param bases
     * @return
     */
    private int getRepeatCount(int n, byte[] bases) {
        if (bases.length < n)
            return 0;
        int count = 1;
        while (true) {
            int pos = count * n;
            if ((pos + n) > bases.length)
                return count;
            for (int i = 0; i < n; i++, pos++) {
                if (bases[i] != bases[pos])
                    return count;
            }
            count++;
        }
    }

    /**
     * Calculates the median absolute deviation for the given statistics.
     *
     * @param statistics
     * @return
     */
    private double getMedianAbsoluteDeviation(DescriptiveStatistics statistics) {
        double median = statistics.getPercentile(50);
        DescriptiveStatistics absoluteDeviations = new DescriptiveStatistics();
        for (int i = 0; i < statistics.getN(); i++) {
            absoluteDeviations.addValue(Math.abs(statistics.getElement(i) - median));
        }
        return absoluteDeviations.getPercentile(50);
    }

    /**
     * Write output VCF file with SNV metrics added.
     *
     * @param metricsLookup the metrics lookup
     */
    private void writeVcf(Map<String, Map<String, Object>> metricsLookup) {
        VCFFileReader reader = new VCFFileReader(inputVcfFile, false);
        VCFHeader header = reader.getFileHeader();
        addInfoHeaderLines(header);

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcfFile)
                .setOutputFileType(OutputType.VCF).setReferenceDictionary(header.getSequenceDictionary())
                .clearOptions();
        VariantContextWriter writer = builder.build();
        writer.writeHeader(header);

        for (VariantContext variant : reader) {
            if (variant.isSNP()) {
                String id = getSnvId(variant);
                Map<String, Object> metrics = metricsLookup.get(id);
                if (metrics == null) {
                    logger.warn("No metrics retrieved from lookup for SNV: " + id);
                } else {
                    variant.getCommonInfo().putAttributes(metrics);
                }
            }
            writer.add(variant);
        }

        writer.close();
        reader.close();
    }

    /**
     * Add header lines for the added INFO fields.
     *
     * @param header the VCF header
     */
    private void addInfoHeaderLines(VCFHeader header) {

        header.addMetaDataLine(new VCFInfoHeaderLine(DEPTH, 1, VCFHeaderLineType.Integer,
                "The number of reads covering the variant position including reads that fall below minimum base and mapping quality thresholds."));
        header.addMetaDataLine(new VCFInfoHeaderLine(DEPTH_CONTROL, 1, VCFHeaderLineType.Integer,
                "The number of reads covering the variant position in the control sample(s) including reads that fall below minimum base and mapping quality thresholds."));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_COUNT, 1, VCFHeaderLineType.Integer,
                "The number of reads covering the variant position excluding reads that fall below minimum base and mapping quality thresholds."));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_COUNT_CONTROL, 1, VCFHeaderLineType.Integer,
                "The number of reads covering the variant position in the control sample(s) excluding reads that fall below minimum base and mapping quality thresholds."));

        header.addMetaDataLine(new VCFInfoHeaderLine(LOW_MAPPING_QUALITY, 1, VCFHeaderLineType.Float,
                "The proportion of all reads from all samples at the variant position that have low mapping quality (less than "
                        + minimumMappingQuality + ")."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_COUNT, 1, VCFHeaderLineType.Integer,
                "The variant allele count, i.e. the number of reads supporting the variant allele."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_COUNT_CONTROL, 1, VCFHeaderLineType.Integer,
                "The variant allele count in the control sample(s)."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_FREQUENCY, 1, VCFHeaderLineType.Float,
                "The variant allele frequency, i.e. the fraction of reads supporting the variant allele."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_FREQUENCY_CONTROL, 1, VCFHeaderLineType.Float,
                "The variant allele frequency in the control sample(s)."));

        header.addMetaDataLine(new VCFInfoHeaderLine(STRAND_BIAS, 1, VCFHeaderLineType.Float,
                "The strand bias for all reads covering the variant position."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_STRAND_BIAS, 1, VCFHeaderLineType.Float,
                "The strand bias for variant-supporting reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(REFERENCE_STRAND_BIAS, 1, VCFHeaderLineType.Float,
                "The strand bias for reference-supporting reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_BASE_QUALITY_MEAN, 1, VCFHeaderLineType.Float,
                "The mean base quality at the variant position of variant reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_BASE_QUALITY_MEDIAN, 1, VCFHeaderLineType.Float,
                "The median base quality at the variant position of variant reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_MAPPING_QUALITY_MEAN, 1, VCFHeaderLineType.Float,
                "The mean mapping quality of variant reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_MAPPING_QUALITY_MEDIAN, 1, VCFHeaderLineType.Float,
                "The median mapping quality of variant reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(MAPPING_QUALITY_MEAN_DIFFERENCE, 1, VCFHeaderLineType.Float,
                "The difference in the mean mapping quality of variant and reference reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(MAPPING_QUALITY_MEDIAN_DIFFERENCE, 1, VCFHeaderLineType.Float,
                "The difference in the median mapping quality of variant and reference reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_MMQS_MEAN, 1, VCFHeaderLineType.Float,
                "The mean mismatch quality sum for variant reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_MMQS_MEDIAN, 1, VCFHeaderLineType.Float,
                "The median mismatch quality sum for variant reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(MMQS_MEAN_DIFFERENCE, 1, VCFHeaderLineType.Float,
                "The difference in mean mismatch quality sum of variant and reference reads."));
        header.addMetaDataLine(new VCFInfoHeaderLine(MMQS_MEDIAN_DIFFERENCE, 1, VCFHeaderLineType.Float,
                "The difference in median mismatch quality sum of variant and reference reads."));

        header.addMetaDataLine(new VCFInfoHeaderLine(DISTANCE_TO_ALIGNMENT_END_MEAN, 1, VCFHeaderLineType.Float,
                "The mean shortest distance of the variant position within the read to either aligned end."));
        header.addMetaDataLine(new VCFInfoHeaderLine(DISTANCE_TO_ALIGNMENT_END_MEDIAN, 1, VCFHeaderLineType.Float,
                "The median shortest distance of the variant position within the read to either aligned end."));
        header.addMetaDataLine(new VCFInfoHeaderLine(DISTANCE_TO_ALIGNMENT_END_MAD, 1, VCFHeaderLineType.Float,
                "The median absolute deviation of the shortest distance of the variant position within the read to either aligned end."));

        header.addMetaDataLine(new VCFInfoHeaderLine(HOMOPOLYMER_LENGTH, 1, VCFHeaderLineType.Integer,
                "The longest continuous homopolymer surrounding or adjacent to the variant position."));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_LENGTH, 1, VCFHeaderLineType.Integer,
                "The length of repetitive sequence adjacent to the variant position where repeats can be 1-, 2-, 3- or 4-mers."));
    }
}
