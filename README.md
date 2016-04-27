# dnaloop
[![Build Status](https://travis-ci.org/aryeelab/dnaloop.svg?branch=master)](https://travis-ci.org/aryeelab/dnaloop)

A preprocessing and QC pipeline for ChIA-PET data.

# Dependencies

The following dependencies need to be manually installed:

- samtools
- bedtools
- bwa

# Installation

Simply run:

    $ pip install dnaloop

# Usage example

The example below uses the test dataset bundled with the `dnaloop` package source code. Download the package with:

`git clone https://github.com/aryeelab/dnaloop.git`


1. Create a tab-separated sample description file with three columns:
  
 - Sample name
 - Read 1 FASTQ
 - Read 2 FASTQ
  
  For example:
  ```
  $ cd dnaloop/test
  $ cat samples.txt 
  naive_esc_1     fastq/naive_esc_1.r1.fastq.gz   fastq/naive_esc_1.r2.fastq.gz
  naive_esc_2     fastq/naive_esc_2.r1.fastq.gz   fastq/naive_esc_2.r2.fastq.gz
  primed_esc_1    fastq/primed_esc_1.r1.fastq.gz  fastq/primed_esc_1.r2.fastq.gz
  primed_esc_2    fastq/primed_esc_2.r1.fastq.gz  fastq/primed_esc_2.r2.fastq.gz
  ```
  
  Note that the FASTQ columns (2 and 3) can contain a comma-separated list of FASTQs. This is commonly the case when a library is sequenced multiple times.
  
2. Run the pipeline:
  ```
  $ preprocess_chiapet --out naive_vs_primed --bwa-index ./test_genome.fa samples.txt
  ```

# Usage details
  ```
  $ preprocess_chiapet --help
  Usage: preprocess_chiapet [OPTIONS] MANIFEST

  A preprocessing and QC pipeline for ChIA-PET data.

  Options:
    --out TEXT         Output directory name  [required]
    --bwa-index TEXT   BWA index location  [required]
    --peak-pad TEXT    Peak padding width (applied on both left and right)
    --merge-gap TEXT   Max gap size for merging peaks
    --use-lsf          Submit jobs to an LSF cluster?
    --bsub-opts TEXT   LSF bsub options
    --keep-temp-files  Keep temporary files?
    --no-qc-report     Skip QC report generation? (Requires R)
    --help             Show this message and exit.

