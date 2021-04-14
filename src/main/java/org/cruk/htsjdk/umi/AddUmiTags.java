/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk.umi;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.cruk.htsjdk.CommandLineProgram;
import org.cruk.htsjdk.ProgressLogger;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Command line tool for adding unique molecular identifier (UMI) tags to
 * records within a BAM file where the UMI or barcode is the first N bases in
 * the read sequence.
 *
 * This utility was written to add UMI tags for sequence data generated from
 * libraries prepared with the Rubicon ThruPLEX Tag-seq kit. The first 6 bases
 * at the beginning of each read constitute the unique molecular identifier in
 * this case. Adding the barcode sequence for each read allows for duplicate
 * marking that takes molecular barcodes into account, e.g. using the Picard
 * MarkDuplicates tool.
 *
 * Note that this tool does not add UMI tags to secondary or supplementary
 * alignment records. These are not normally marked as duplicates by Picard
 * MarkDuplicates where the data are in coordinate sort order.
 *
 * @author eldrid01
 */
public class AddUmiTags extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger();

    private File inputBamFile;
    private File outputBamFile;
    private int umiLength;
    private String umiTag;

    public static void main(String[] args) {
        AddUmiTags addUmiTags = new AddUmiTags();
        addUmiTags.parseCommandLineArgs(args);
        addUmiTags.run();
    }

    @Override
    protected String getHelpDescription() {
        return "Adds SAM tags for UMI sequences found in the specified number of bases at the beginning of each read.";
    }

    @Override
    protected Options createOptions() {
        Options options = super.createOptions();

        Option option = new Option("i", "input", true, "Input BAM file (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("o", "output", true, "Output BAM file (required)");
        option.setRequired(true);
        option.setType(File.class);
        options.addOption(option);

        option = new Option("l", "umi-length", true,
                "The number of bases at the beginning of the read that consitute the UMI or barcode (required)");
        option.setRequired(true);
        option.setType(Number.class);
        options.addOption(option);

        option = new Option("t", "umi-tag", true,
                "The tag to use for the UMI or barcode in the output BAM file (required)");
        option.setRequired(true);
        options.addOption(option);

        return options;
    }

    @Override
    protected void extractOptionValues(CommandLine commandLine) throws ParseException {
        inputBamFile = (File) commandLine.getParsedOptionValue("input");
        outputBamFile = (File) commandLine.getParsedOptionValue("output");
        umiLength = ((Number) commandLine.getParsedOptionValue("umi-length")).intValue();
        umiTag = commandLine.getOptionValue("umi-tag");
    }

    /**
     * Main run method for iterating over records in a BAM file and adding a tag for
     * the unique molecular index (UMI) tag found in the first specified number of
     * bases.
     */
    private void run() {
        ProgressLogger progress = new ProgressLogger(logger, 1000000);

        IOUtil.assertFileIsReadable(inputBamFile);
        IOUtil.assertFileIsWritable(outputBamFile);

        SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
                .open(inputBamFile);

        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeSAMOrBAMWriter(reader.getFileHeader(), true, outputBamFile);

        SAMRecordIterator iterator = reader.iterator();

        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();

            if (!record.isSecondaryOrSupplementary()) {
                String sequence = record.getReadString();
                if (record.getReadNegativeStrandFlag()) {
                    sequence = SequenceUtil.reverseComplement(sequence);
                }
                String umi = sequence.substring(0, umiLength);
                record.setAttribute(umiTag, umi);
            }

            writer.addAlignment(record);

            progress.record(record);
        }

        iterator.close();

        logger.info("Writing " + outputBamFile.getName());
        writer.close();

        CloserUtil.close(reader);
    }
}
