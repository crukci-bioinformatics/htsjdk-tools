# htsjdk-tools

Command-line tools for processing high-throughput sequencing (HTS) data and
genomic variants using the HTSJDK library.

This package contains some command-line utilities for processing high-throughput
sequencing (HTS) data and genomic variants using the
[HTSJDK](https://github.com/samtools/htsjdk) Java library.

_htsjdk-tools_ contains the following utilities:

* **pileup-counts** - Generate a pileup summary with read counts for each position and allele
* **calculate-snv-metrics** - Calculate metrics for single nucleotide variants (SNVs) from a VCF file based on aligned reads within the given BAM file(s)
* **add-umi-tags** - Add SAM tags for UMI sequences found in the specified number of bases at the beginning of each read

## Installation

### Pre-requisites

_htsjdk-tools_ requires a Java runtime, e.g. Java SE 8 or above.

### Installing a pre-packaged release

Packaged releases of _htsjdk-tools_ can be downloaded from the
[releases](https://github.com/crukci-bioinformatics/htsjdk-tools/releases) page.

Unpack the downloaded tarball to an installation directory, substituting
the version number as appropriate.

        tar zxf htsjdk-tools-1.0-distribution.tar.gz

The command-line utilities are located in the _bin_ subdirectory
(htsjdk-tools-1.0/bin).

### Building from source

The _htsjdk-tools_ package is built using Apache Maven, a software project
management and build automation tool. Details on how to install and run Maven
can be found [here](http://maven.apache.org).

1. Clone the project

        git clone https://github.com/crukci-bioinformatics/htsjdk-tools.git
        cd htsjdk-tools

2. Build the package

        mvn package

3. Unpack the _htsjdk-tools_ tarball to an installation directory, substituting
the version number as appropriate

        tar zxf target/htsjdk-tools-1.0-SNAPSHOT-distribution.tar.gz

This will create a directory named htsjdk-tools-1.0-SNAPSHOT which can be moved
to a preferred installation location.

## Running _htsjdk-tools_

Details about how to run each of the tools and the command-line options or
parameter settings that are available are listed by running the tool with the
`--help` argument, e.g.

        pileup-counts --help

## SNV metrics

The `calculate-snv-metrics` utility computes various metrics for single
nucleotide variants (SNVs) from a given VCF file based on the variant-supporting
and reference sequence reads within the specified input BAM file(s).

Values for the metrics computed are added to the INFO column of the output VCF
file.

Usage example:

        calculate-snv-metrics \
          --reference-sequence homo_sapiens.fa \
          --variants snv.vcf \
          --input my.bam
          --output snv.metrics.vcf

This will annotate the SNV variants within the input VCF file, `snv.vcf`, based
on reads within the input BAM file. Each of the various SNV metrics are
described below.

Some metrics are computed from variant-supporting reads and others for reads in
which the base aligned matches the reference sequence. The tool allows for
filtering of reads used in computing these metrics based on the sample
identifier in the SM tag within the SAM record. This can be useful when
calculating metrics for somatic variant calls from cancer genome sequencing
where both the tumour sample and a matched normal were sequenced. In this
scenario, the metrics for variant-supporting reads can be restricted to just the
sequence data from the tumour sample while the metrics for reference sequence
reads can be calculated solely for data from the matched normal sample by
specifying the SM identifiers for each using `--sample` and `--control-sample`
options.

        calculate-snv-metrics \
          --reference-sequence homo_sapiens.fa \
          --variants snv.vcf \
          --input tumour.bam \
          --input control.bam \
          --sample tumour123 \
          --control-sample control123 \
          --minimum-mapping-quality 1 \
          --minimum-base-quality 10 \
          --output snv.metrics.vcf

Two input BAM files are provided in this scenario for somatic variant calling,
one for the tumour sample and the other for the matched normal (control) sample,
but `calculate-snv-metrics` using the SM tags to determine which reads to use
for each metric. Multiple sample identifiers can be specified if necessary.

The following metrics are computed for each SNV:

Metric | Description
-------|------------
Depth |The number of reads covering the variant position including duplicates and reads that fall below minimum base and mapping quality thresholds.
ReadCount | The number of reads covering the variant position excluding duplicates and reads that fall below minimum base and mapping quality thresholds.
VariantAlleleCount | The variant allele count, i.e. the number of reads supporting the variant allele.
VariantAlleleFrequency | The variant allele frequency, i.e. the fraction of reads supporting the variant allele.
DepthControl | The number of reads covering the variant position in the control sample(s) including duplicates and reads that fall below minimum base and mapping quality thresholds
ReadCountControl | The number of reads covering the variant position in the control sample(s) excluding duplicates and reads that fall below minimum base and mapping quality thresholds.
VariantAlleleCountControl | The variant allele count in the control sample(s).
StrandBias | The strand bias for all reads covering the variant position.
VariantStrandBias | The strand bias for variant-supporting reads.
ReferenceStrandBias | The strand bias for reference-supporting reads.
LowMapQual | The proportion of all reads from all samples at the variant position that have low mapping quality (less than the specified minimumMappingQuality).
VariantBaseQual | The mean base quality at the variant position of variant reads.
VariantBaseQualMedian | The median base quality at the variant position of variant reads.
VariantMapQual | The mean mapping quality of variant reads.
VariantMapQualMedian | The median mapping quality of variant reads.
MapQualDiff | The difference in the mean mapping quality of variant and reference reads.
MapQualDiffMedian | The difference in the median mapping quality of variant and reference reads.
VariantMMQS | The mean mismatch quality sum for variant reads.
VariantMMQSMedian | The median mismatch quality sum for variant reads.
MMQSDiff | The difference in mean mismatch quality sum of variant and reference reads.
MMQSDiffMedian | The difference in median mismatch quality sum of variant and reference reads.
DistanceToAlignmentEnd | The mean shortest distance of the variant position within the read to either aligned end.
DistanceToAlignmentEndMedian | The median shortest distance of the variant position within the read to either aligned end.
DistanceToAlignmentEndMAD | The median absolute deviation of the shortest distance of the variant position within the read to either aligned end.
HomopolymerLength | The longest continuous homopolymer surrounding or adjacent to the variant position.
Repeat | The length of repetitive sequence adjacent to the variant position where repeats can be 1-, 2-, 3- or 4-mers.

### Filtering based on SNV metrics

The VariantFiltration tool from [GATK](https://gatk.broadinstitute.org) can be
used to filter variants based on the SNV metrics computed using
`calculate-snv-metrics`.

Based on alignment with BWA MEM for whole genome sequencing datasets, the
following filters are recommended for MuTect2 and Strelka variant callers.

These filters have been tuned using the ICGC benchmark datasets from Alioto et
al., Nature Commun. 2015, 6:10001 (http://www.ncbi.nlm.nih.gov/pubmed/26647970)
and tested on synthetic datasets from the ICGC-TCGA DREAM Mutation Calling
challenge (https://www.synapse.org/#!Synapse:syn312572).

MuTect2 filters (more stringent for higher precision):

    gatk VariantFiltration \
      --variant snv.metrics.vcf \
      --output snv.filtered.vcf \
      --filter-name VariantAlleleCount \
      --filter-expression "VariantAlleleCount < 4" \
      --filter-name VariantAlleleCountControl \
      --filter-expression "VariantAlleleCountControl > 1" \
      --filter-name VariantMapQualMedian \
      --filter-expression "VariantMapQualMedian < 40.0" \
      --filter-name MapQualDiffMedian \
      --filter-expression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
      --filter-name LowMapQual \
      --filter-expression "LowMapQual > 0.05" \
      --filter-name VariantBaseQualMedian \
      --filter-expression "VariantBaseQualMedian < 30.0" \
      --filter-name StrandBias \
      --filter-expression "VariantAlleleCount >= 7 && VariantStrandBias < 0.05 && ReferenceStrandBias >= 0.2"

MuTect2 filters (less stringent for higher recall at expense of precision):

    gatk VariantFiltration \
      --variant snv.metrics.vcf \
      --output snv.filtered.vcf \
      --filter-name VariantAlleleCount \
      --filter-expression "VariantAlleleCount < 3" \
      --filter-name VariantAlleleCountControl \
      --filter-expression "VariantAlleleCountControl > 1" \
      --filter-name VariantMapQualMedian \
      --filter-expression "VariantMapQualMedian < 40.0" \
      --filter-name MapQualDiffMedian \
      --filter-expression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
      --filter-name LowMapQual \
      --filter-expression "LowMapQual > 0.05" \
      --filter-name VariantBaseQualMedian \
      --filter-expression "VariantBaseQualMedian < 25.0"

Strelka filters (includes position in read filter):

    gatk VariantFiltration \
      --variant snv.metrics.vcf \
      --output snv.filtered.vcf \
      --filter-name VariantAlleleCount \
      --filter-expression "VariantAlleleCount < 4" \
      --filter-name VariantAlleleCountControl \
      --filter-expression "VariantAlleleCountControl > 1" \
      --filter-name VariantMapQualMedian \
      --filter-expression "VariantMapQualMedian < 40.0" \
      --filter-name MapQualDiffMedian \
      --filter-expression "MapQualDiffMedian < -5.0 || MapQualDiffMedian > 5.0" \
      --filter-name LowMapQual \
      --filter-expression "LowMapQual > 0.05" \
      --filter-name VariantBaseQualMedian \
      --filter-expression "VariantBaseQualMedian < 30.0" \
      --filter-name StrandBias \
      --filter-expression "VariantAlleleCount >= 7 && VariantStrandBias < 0.05 && ReferenceStrandBias >= 0.2" \
      --filter-name DistanceToAlignmentEndMedian \
      --filter-expression "DistanceToAlignmentEndMedian < 10.0" \
      --filter-name DistanceToAlignmentEndMAD \
      --filter-expression "DistanceToAlignmentEndMAD < 3.0"
