#!/bin/bash

#source /apps/lab/aryee/pyenv/versions/venv-2.7.10/bin/activate

# Argument parsing can be fixed later, either by making this a Python script,
# or through more elegant bash. (e.g. http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash)
SAMPLE_DIR=$1
BWA_INDEX=$2
MERGE_GAP=$3
R1_FASTQ=$4
R2_FASTQ=$5

## Test argument set 1
#SAMPLE_DIR="output/tiny"
#R1_FASTQ="input/test/tiny_r1.fastq.gz"
#R2_FASTQ="input/test/tiny_r1.fastq.gz"

## Test argument set 2. This combines two input fastqs for R1, R2
#SAMPLE_DIR="small_and_tiny"
#R1_FASTQ="small_r1.fastq.gz,tiny_r1.fastq.gz"
#R2_FASTQ="small_r2.fastq.gz,tiny_r1.fastq.gz"

# Parameters are hard-coded for now:
MIN_QUAL=30
READ_LEN=100 # Only used to define the intervals in BEDPE (i.e. start=map_pos, end=map_pos+READ_LEN)
FWD_ADAPTER="ACGCGATATCTTATCTGACT"
REV_ADAPTER="AGTCAGATAAGATATCGCGT"
#BWA_INDEX="/data/aryee/pub/genomes/grch37/bwa_index/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
#BWA_INDEX="`pwd`/test/test_genome.fa"

mkdir -p $SAMPLE_DIR

LOG_FILE="preprocess_fastq.log"
echo "`date`: Run starting in $SAMPLE_DIR" | tee $SAMPLE_DIR/$LOG_FILE
echo "`date`: Copying ${R1_FASTQ} to r1.fastq.gz" | tee -a $SAMPLE_DIR/$LOG_FILE
# Note:  ${VAR//,/ } replaces commas in VAR with spaces
cat ${R1_FASTQ//,/ } > $SAMPLE_DIR/r1.fastq.gz
echo "`date`: Copying ${R2_FASTQ} to r2.fastq.gz" | tee -a $SAMPLE_DIR/$LOG_FILE
cat ${R2_FASTQ//,/ } > $SAMPLE_DIR/r2.fastq.gz
cd $SAMPLE_DIR

echo "`date`: Finding and removing linkers" | tee -a $LOG_FILE
cutadapt -n 3 -m 17 --overlap 10 --pair-filter=both --suffix " {name}" -a forward=$FWD_ADAPTER -a reverse=$REV_ADAPTER -A forward=$FWD_ADAPTER -A reverse=$REV_ADAPTER -o r1.linker_removed.fastq.gz -p r2.linker_removed.fastq.gz --untrimmed-output r1.no_linker.fastq.gz --untrimmed-paired-output r2.no_linker.fastq.gz r1.fastq.gz r2.fastq.gz > cutadapt.log

echo "`date`: Writing linker counts to linker_stats.txt and linker_stats_detail.txt" | tee -a $LOG_FILE
zcat r1.linker_removed.fastq.gz | awk 'NR%4==1' | awk '{print $NF}' > linker_r1.txt
zcat r2.linker_removed.fastq.gz | awk 'NR%4==1' | awk '{print $NF}' > linker_r2.txt
paste linker_r1.txt linker_r2.txt > linker.txt
cat linker.txt | sort | uniq -c | sed 's/adapter/linker/' > linker_stats_detail.txt
sed 's/forward/linker/g; s/reverse/linker/g; s/adapter/linker/' linker.txt | sort | uniq -c > linker_stats.txt 
rm linker.txt
rm linker_r1.txt
rm linker_r2.txt

echo "`date`: Aligning PETs that had a linker (R1 or R2 or both)" | tee -a $LOG_FILE
bwa mem $BWA_INDEX r1.linker_removed.fastq.gz 2> bwa.log  | samtools view -bS - > r1.bam 
bwa mem $BWA_INDEX r2.linker_removed.fastq.gz 2>> bwa.log | samtools view -bS - > r2.bam 

echo "`date`: Writing interactions with both anchors mapped (MAPQ>=$MIN_QUAL) to interactions.bedpe" | tee -a $LOG_FILE
samtools view -F2304 r1.bam  | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > pos_r1.bed
samtools view -F2304 r2.bam  | awk -v READ_LEN="$READ_LEN" -v MIN_QUAL="$MIN_QUAL" '{if ($5>=MIN_QUAL) print $3,$4,$4+READ_LEN,$1; else print "*","*","*",$1}' OFS='\t' > pos_r2.bed
# Make sure reads match up between R1 and R2
cut -f 4 pos_r1.bed > names_r1.txt
cut -f 4 pos_r2.bed > names_r2.txt
DIFF=$(diff names_r1.txt names_r2.txt) 
if [ "$DIFF" != "" ] 
then
    echo "ERROR: Read names don't match between pos_r1.bed and pos_r2.bed"
    exit 1
fi

paste pos_r1.bed pos_r2.bed | cut -f 1-3,5-8 > interactions.tmp
# Keep only interactions where both reads are mapped
awk '{if ($1 != "*" && $4 != "*") print}' interactions.tmp > interactions.unsorted.bedpe
# Swap anchors if necessary to make the position of anchor1 < anchor2
awk '{if ($1<$4 || ($1==$4 && $2<=$5)) print $0; else print $4,$5,$6,$1,$2,$3,$7}' 'OFS=\t' interactions.unsorted.bedpe > interactions.unsorted.orderedanchors.bedpe
sort -k1,1V -k2,2n -k4,4V -k5,5n interactions.unsorted.orderedanchors.bedpe > interactions.bedpe
echo "`date`: Writing unique (i.e. deduplicated) interactions to interactions.dedup.bedpe (based on chr1, pos1, chr2, pos2)" | tee -a $LOG_FILE
sort -k1,1V -k2,2n -k4,4V -k5,5n --unique interactions.bedpe > interactions.dedup.bedpe

echo "`date`: Writing read processing summary counts to read_stats.txt" | tee -a $LOG_FILE
STATS_FILE='read_stats.txt'
echo "Sample_directory: $SAMPLE_DIR" > $STATS_FILE
FASTQ_LEN=`zcat r1.fastq.gz | wc -l`
echo "Total_PETs: $((FASTQ_LEN / 4))" >> $STATS_FILE
LINKER_FASTQ_LEN=`zcat r1.linker_removed.fastq.gz | wc -l`
echo "PETs_with_linker: $((LINKER_FASTQ_LEN / 4))" >> $STATS_FILE
NO_LINKER_FASTQ_LEN=`zcat r1.no_linker.fastq.gz | wc -l`
echo "PETs_without_linker: $((NO_LINKER_FASTQ_LEN / 4))" >> $STATS_FILE
echo "Discarded_short_PETs: $(((FASTQ_LEN - LINKER_FASTQ_LEN - NO_LINKER_FASTQ_LEN) / 4))" >> $STATS_FILE
echo "Linker_on_R1_and_R2: `cat linker_stats.txt | grep -w $'linker\tlinker' | awk '{print $1}'`" >> $STATS_FILE
echo "Linker_on_R1_only: `cat linker_stats.txt | grep -w $'linker\tno_linker' | awk '{print $1}'`" >> $STATS_FILE
echo "Linker_on_R2_only: `cat linker_stats.txt | grep -w $'no_linker\tlinker' | awk '{print $1}'`" >> $STATS_FILE
echo "Mapped_r1_q5: `samtools view -q5 -F2304 r1.bam | wc -l`" >> $STATS_FILE
echo "Mapped_r2_q5: `samtools view -q5 -F2304 r2.bam | wc -l`" >> $STATS_FILE
echo "Mapped_r1_q${MIN_QUAL}: `samtools view -q${MIN_QUAL} -F2304 r1.bam | wc -l`" >> $STATS_FILE
echo "Mapped_r2_q${MIN_QUAL}: `samtools view -q${MIN_QUAL} -F2304 r2.bam | wc -l`" >> $STATS_FILE
echo "Mapped_PETs_q${MIN_QUAL}: `cat interactions.bedpe | wc -l`" >> $STATS_FILE
echo "Mapped_unique_PETs_q${MIN_QUAL}: `cat interactions.dedup.bedpe | wc -l`" >> $STATS_FILE
echo "Mapped_intrachromosal_pairs_q${MIN_QUAL}: `awk '{if ($1==$4) print}' interactions.bedpe | wc -l`" >> $STATS_FILE
echo "Mapped_unique_intrachromosal_pairs_q${MIN_QUAL}: `awk '{if ($1==$4) print}' interactions.dedup.bedpe | wc -l`" >> $STATS_FILE
echo "Mapped_unique_intrachromosal_pairs_q${MIN_QUAL}_5kb: `awk '{if (($1==$4) && ($5-$3 >= 5000)) print}' interactions.dedup.bedpe | wc -l`" >> $STATS_FILE


##############################################################################
####  Note that the section below dealing with peak calling and
####  anchor mapping largely duplicates code from preprocess_chiapet_set. DRY
echo "`date`: Running preprocess_chiapet_set to compute anchor stats" | tee -a $LOG_FILE
SAMPLE_ANCHOR_DIR=anchors
mkdir -p $SAMPLE_ANCHOR_DIR
echo "`date`: Writing unique interaction left and right anchors for $SAMPLE_DIR to left.dedup.bed and right.dedup.bed" | tee -a $LOG_FILE
cut -f 1-3,7 interactions.dedup.bedpe > $SAMPLE_ANCHOR_DIR/left.dedup.bed 
cut -f 4-6,7 interactions.dedup.bedpe > $SAMPLE_ANCHOR_DIR/right.dedup.bed 
cd $SAMPLE_ANCHOR_DIR
echo "`date`: Adding reads to reads.bed" | tee -a ../$LOG_FILE
cat left.dedup.bed right.dedup.bed > reads.bed
echo "`date`: Calling peaks using reads.bed to define interaction anchor locations" | tee -a ../$LOG_FILE
echo "`date`: Using `which macs2`" | tee -a ../$LOG_FILE
echo "`date`: MACS command: macs2 callpeak -t reads.bed -f BED -n anchor --nomodel -p 0.01 --outdir ." | tee -a ../$LOG_FILE
macs2 callpeak -t reads.bed -f BED -n anchor --nomodel -p 0.01 --outdir .
NUM_PEAKS=`cat anchor_peaks.narrowPeak | wc -l`
echo "`date`: Found $NUM_PEAKS peaks" | tee -a ../$LOG_FILE
bedtools merge -d $MERGE_GAP -i anchor_peaks.narrowPeak > anchor_peaks.merged.bed
NUM_MERGED_PEAKS=`cat anchor_peaks.merged.bed | wc -l`
echo "`date`: Merged peaks within ${MERGE_GAP}bp resulting in $NUM_MERGED_PEAKS anchors" | tee -a ../$LOG_FILE
echo "`date`:   Mapping `cat left.dedup.bed | wc -l` PETs to anchors by intersecting with peaks" | tee -a ../$LOG_FILE
bedtools intersect -loj -a left.dedup.bed -b anchor_peaks.merged.bed | awk '{print $5,$6,$7,$4}' OFS='\t' > anchor1.bed
bedtools intersect -loj -a right.dedup.bed -b anchor_peaks.merged.bed | awk '{print $5,$6,$7,$4}' OFS='\t' > anchor2.bed
# Confirm reads match up between left and right
cut -f 4 anchor1.bed > anchor1_names.txt
cut -f 4 anchor2.bed > anchor2_names.txt
DIFF=$(diff anchor1_names.txt anchor1_names.txt) 
if [ "$DIFF" != "" ] 
then
    echo "ERROR: Read names don't match between anchor1.bed and anchor2.bed" | tee -a ../$LOG_FILE
    exit 1
fi
# Create BEDPE
paste anchor1.bed anchor2.bed | cut -f 1-3,5-8 > anchor_interactions.tmp
# Keep only interactions where both reads map to anchors
awk '{if ($1 != "." && $4 != ".") print}' anchor_interactions.tmp > anchor_interactions.bedpe
ANCHOR_PETS=`cat anchor_interactions.bedpe | wc -l`
echo "`date`:   Wrote $ANCHOR_PETS PETs where both reads map to anchors (anchor_interactions.bedpe)" | tee -a ../$LOG_FILE    
DIFF_ANCHOR_PETS_5KB=`awk '($1==$4) && ($5-$3 >= 5000)' anchor_interactions.tmp | wc -l`    
cut -f1-6 anchor_interactions.bedpe | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' > loop_counts.bedpe
SAME_ANCHOR_LOOPS=`cat loop_counts.bedpe | awk '$1==$4 && $2==$5' | wc -l`
DIFF_ANCHOR_LOOPS=`cat loop_counts.bedpe | awk '$1!=$4 || $2!=$5' | wc -l`
DIFF_ANCHOR_LOOPS_3PETS=`cat loop_counts.bedpe | awk '($1!=$4 || $2!=$5) && $8>=3' | wc -l`
DIFF_ANCHOR_INTRACHROMOSOMAL_LOOPS_3PETS=`cat loop_counts.bedpe | awk '($1==$4) && ($2!=$5) && $8>=3' | wc -l`
DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_1PET=`cat loop_counts.bedpe | awk '($1==$4) && ($5-$3 >= 5000)' | wc -l`
DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_2PETS=`cat loop_counts.bedpe | awk '($1==$4) && ($5-$3 >= 5000) && $8>=2' | wc -l`
DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_3PETS=`cat loop_counts.bedpe | awk '($1==$4) && ($5-$3 >= 5000) && $8>=3' | wc -l`
echo "Num_peaks: $NUM_PEAKS" >> ../$STATS_FILE
echo "Num_merged_peaks: $NUM_MERGED_PEAKS" >> ../$STATS_FILE
echo "Anchor_mapped_pets: $ANCHOR_PETS" >> ../$STATS_FILE
echo "Anchor_mapped_pets_5kb: $DIFF_ANCHOR_PETS_5KB" >> ../$STATS_FILE
echo "Same_anchor_loops: $SAME_ANCHOR_LOOPS" >> ../$STATS_FILE
echo "5kb_loops_1pet: $DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_1PET" >> ../$STATS_FILE
echo "5kb_loops_2pets: $DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_2PETS" >> ../$STATS_FILE
echo "5kb_loops_3pets: $DIFF_ANCHOR_INTRACHROMOSOMAL_5KB_LOOPS_3PETS" >> ../$STATS_FILE

cd ..
echo "`date`: Run finished. See $STATS_FILE for statistics" | tee -a $LOG_FILE
