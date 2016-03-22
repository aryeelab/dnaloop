# dnaloop

A preprocessing and QC pipeline for ChIA-PET data.


# Installation

Use [`pipsi`](https://github.com/mitsuhiko/pipsi#readme) to install.

Simply run:

    $ pipsi install .


# Usage example

**Step 1: Create a tab-separated sample description file with three columns**
  
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
  
**Step 2: Run the pipeline**

    $ cd test
    $ preprocess_chiapet --out naive_vs_primed --bwa-index test_genome.fa samples.txt

**Output**

The key output files are:

- `[SAMPLE_NAME].loop_counts.bedpe`
- `preprocess_chiapet_set.log`

e.g.

    $ head naive_vs_primed/naive_esc_1.loop_counts.bedpe 
    chrZ 19843 20399 chrZ 19843 20399 . 14
    chrZ 19843 20399 chrZ 21941 22516 . 1
    chrZ 19843 20399 chrZ 40625 41746 . 3
    chrZ 21941 22516 chrZ 21941 22516 . 10
    chrZ 21941 22516 chrZ 40625 41746 . 1


# Usage details

    $ preprocess_chiapet --help
    Usage: preprocess_chiapet [OPTIONS] MANIFEST

      A preprocessing and QC pipeline for ChIA-PET data.

    Options:
      --out TEXT        Output folder name  [required]
      --bwa-index TEXT  BWA Index  [required]
      --help            Show this message and exit.

