# dnaloop

A preprocessing and QC pipeline for ChIA-PET data.


# Installation

Use [`pipsi`](https://github.com/mitsuhiko/pipsi#readme) to install.

Simply run:

    $ pipsi install .


# Usage example

1. Create a sample description manifest with three columns:

- Sample name
- Read 1 FASTQ 
- Read 2 FASTQ

For example:

    $ cat test/samples.txt 
    naive_esc_1     fastq/naive_esc_1.r1.fastq.gz   fastq/naive_esc_1.r2.fastq.gz
    naive_esc_2     fastq/naive_esc_2.r1.fastq.gz   fastq/naive_esc_2.r2.fastq.gz
    primed_esc_1    fastq/primed_esc_1.r1.fastq.gz  fastq/primed_esc_1.r2.fastq.gz
    primed_esc_2    fastq/primed_esc_2.r1.fastq.gz  fastq/primed_esc_2.r2.fastq.gz

Note that the FASTQ columns (2 and 3) can contain a comma-separated list of FASTQs. This is common when a sample is sequenced multiple times.

2. Run the pipeline:

    $ cd test
    $ preprocess_chiapet --out naive_vs_primed --bwa-index test_genome.fa samples.txt


