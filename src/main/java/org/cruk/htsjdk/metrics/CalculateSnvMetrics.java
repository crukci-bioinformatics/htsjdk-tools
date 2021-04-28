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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;
import org.cruk.htsjdk.pileup.PileupUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
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

    @Option(names = { "-i",
            "--input" }, required = true, description = "Input BAM file(s) which must be in coordinate sort order and indexed (required).")
    private List<File> bamFiles;

    @Option(names = { "-v",
            "--variants" }, required = true, description = "Input VCF file containing single nucleotide variants (required).")
    private File inputVcfFile;

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

        // metrics lookup
        Map<String, Map<String, Object>> metricsLookup = new HashMap<>();

        while (mergingSamLocusIterator.hasNext()) {
            LocusInfo locusInfo = mergingSamLocusIterator.next();

            List<RecordAndOffset> pileup = locusInfo.getRecordAndOffsets();

            boolean foundMatchingSnv = false;

            while (isVariantAtLocus(snvIterator.peek(), locusInfo)) {
                foundMatchingSnv = true;

                VariantContext snv = snvIterator.next();

                String key = getSnvKey(snv);

                // shouldn't really have multiple entries for the same variant
                // but only compute the metrics once for each distinct SNV
                if (!metricsLookup.containsKey(key)) {
                    Map<String, Object> metrics = calculateSnvMetrics(snv, pileup);
                    metricsLookup.put(key, metrics);
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
                    logger.warn("Warning multiple alternate alleles for variant " + getSnvKey(variant)
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
    private String getSnvKey(VariantContext snv) {
        return String.format("%s:%d %c>%c", snv.getContig(), snv.getStart(), snv.getReference().getBases()[0],
                snv.getAltAlleleWithHighestAlleleCount().getBases()[0]);
    }

    /**
     * Calculates SNV metrics for the given variant and read pileup.
     *
     * For multi-allelic variants, metrics are only calculated for the alternate
     * allele with the highest allele count.
     *
     * @param snv    the SNV variant
     * @param pileup the read pileup
     * @return the computed SNV metrics
     */
    private Map<String, Object> calculateSnvMetrics(VariantContext snv, List<RecordAndOffset> pileup) {

        Map<String, Object> metrics = new HashMap<>();

        byte referenceBase = snv.getReference().getBases()[0];
        byte variantAllele = snv.getAltAlleleWithHighestAlleleCount().getBases()[0];

        // depth metrics calculated before applying additional mapping quality and
        // overlap filters
        metrics.put(DEPTH, samples.isEmpty() ? pileup.size() : PileupUtils.filterSamples(pileup, samples).size());

        if (!controlSamples.isEmpty()) {
            metrics.put(DEPTH_CONTROL, PileupUtils.filterSamples(pileup, controlSamples).size());
        }

        int unfilteredReadCount = pileup.size();

        // apply mapping quality filter
        pileup = PileupUtils.filterLowMappingQualities(pileup, minimumMappingQuality);

        if (unfilteredReadCount > 0) {
            metrics.put(LOW_MAPPING_QUALITY, 1.0 - (pileup.size() / (double) unfilteredReadCount));
        }

        // apply base quality filter
        pileup = PileupUtils.filterLowBaseQualities(pileup, minimumBaseQuality);

        // apply filter for overlapping fragments
        pileup = PileupUtils.filterOverlaps(pileup);

        // apply filter for the specified sample(s)
        List<RecordAndOffset> samplePileup = samples.isEmpty() ? pileup : PileupUtils.filterSamples(pileup, samples);

        // collect reference metrics from control samples if these exist,
        // otherwise from reference-supporting reads in the same set of
        // samples used for computing the variant metrics
        List<RecordAndOffset> referencePileup = controlSamples.isEmpty() ? samplePileup
                : PileupUtils.filterSamples(pileup, controlSamples);

        int readCount = samplePileup.size();
        metrics.put(READ_COUNT, readCount);

        int referenceReadCount = referencePileup.size();
        if (!controlSamples.isEmpty()) {
            metrics.put(READ_COUNT_CONTROL, referenceReadCount);
        }

        int referenceBaseIndex = PileupUtils.getBaseIndex(referenceBase);
        if (referenceBaseIndex == -1) {
            logger.warn("Unrecognized reference base for SNV: " + getSnvKey(snv));
            return metrics;
        }

        int variantAlleleIndex = PileupUtils.getBaseIndex(variantAllele);
        if (variantAlleleIndex == -1) {
            logger.warn("Unrecognized alternate allele for SNV: " + getSnvKey(snv));
            return metrics;
        }

        int[] sampleBaseCounts = PileupUtils.getBaseCounts(samplePileup);
        int[] referenceBaseCounts = PileupUtils.getBaseCounts(referencePileup);

        int variantAlleleCount = sampleBaseCounts[variantAlleleIndex];
        metrics.put(VARIANT_ALLELE_COUNT, variantAlleleCount);

        if (readCount > 0) {
            metrics.put(VARIANT_ALLELE_FREQUENCY, variantAlleleCount / (double) readCount);
        }

        if (!controlSamples.isEmpty()) {
            int variantAlleleCountControl = referenceBaseCounts[variantAlleleIndex];
            metrics.put(VARIANT_ALLELE_COUNT_CONTROL, variantAlleleCountControl);
            if (readCount > 0) {
                metrics.put(VARIANT_ALLELE_FREQUENCY_CONTROL, variantAlleleCountControl / (double) referenceReadCount);
            }
        }

        return metrics;
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
                String key = getSnvKey(variant);
                Map<String, Object> metrics = metricsLookup.get(key);
                if (metrics == null) {
                    logger.warn("No metrics retrieved from lookup for SNV: " + key);
                } else {
                    variant.getCommonInfo().putAttributes(metrics);
                }
            }
            writer.add(variant);
        }

        writer.close();
        reader.close();
    }

    private static final String DEPTH = "Depth";
    private static final String DEPTH_CONTROL = "DepthControl";
    private static final String LOW_MAPPING_QUALITY = "LowMapQual";
    private static final String READ_COUNT = "ReadCount";
    private static final String READ_COUNT_CONTROL = "ReadCountControl";
    private static final String VARIANT_ALLELE_COUNT = "VariantAlleleCount";
    private static final String VARIANT_ALLELE_COUNT_CONTROL = "VariantAlleleCountControl";
    private static final String VARIANT_ALLELE_FREQUENCY = "VariantAlleleFrequency";
    private static final String VARIANT_ALLELE_FREQUENCY_CONTROL = "VariantAlleleFrequencyControl";

    /**
     * Add header lines for the added INFO fields.
     *
     * @param header the VCF header
     */
    private void addInfoHeaderLines(VCFHeader header) {

        header.addMetaDataLine(new VCFInfoHeaderLine(DEPTH, 1, VCFHeaderLineType.String,
                "The number of reads covering the variant position including reads that fall below minimum base and mapping quality thresholds."));
        header.addMetaDataLine(new VCFInfoHeaderLine(DEPTH_CONTROL, 1, VCFHeaderLineType.String,
                "The number of reads covering the variant position in the control sample(s) including reads that fall below minimum base and mapping quality thresholds."));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_COUNT, 1, VCFHeaderLineType.String,
                "The number of reads covering the variant position excluding reads that fall below minimum base and mapping quality thresholds."));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_COUNT_CONTROL, 1, VCFHeaderLineType.String,
                "The number of reads covering the variant position in the control sample(s) excluding reads that fall below minimum base and mapping quality thresholds."));

        header.addMetaDataLine(new VCFInfoHeaderLine(LOW_MAPPING_QUALITY, 1, VCFHeaderLineType.String,
                "The proportion of all reads from all samples at the variant position that have low mapping quality (less than "
                        + minimumMappingQuality + ")."));

        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_COUNT, 1, VCFHeaderLineType.String,
                "The variant allele count, i.e. the number of reads supporting the variant allele."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_COUNT_CONTROL, 1, VCFHeaderLineType.String,
                "The variant allele count in the control sample(s)."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_FREQUENCY, 1, VCFHeaderLineType.String,
                "The variant allele frequency, i.e. the fraction of reads supporting the variant allele."));
        header.addMetaDataLine(new VCFInfoHeaderLine(VARIANT_ALLELE_FREQUENCY_CONTROL, 1, VCFHeaderLineType.String,
                "The variant allele frequency in the control sample(s)."));

        // header.addMetaDataLine(new VCFInfoHeaderLine(ANNOTATION, 1,
        // VCFHeaderLineType.String, ""));
    }
}