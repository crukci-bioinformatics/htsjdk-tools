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
