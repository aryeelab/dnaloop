#!/bin/bash

# Argument parsing can be fixed later, either by making this a Python script,
# or through more elegant bash. (e.g. http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash)
CHIAPET_SET_DIR=$1
MERGE_GAP=$2
SAMPLE_DIRS=${@:3}

#CHIAPET_SET_DIR="output/chiapet_set_1"
#SAMPLE_DIRS="output/tiny output/small"
#SAMPLE_DIRS="output/naive_esc_1 output/primed_esc_1"

echo "CHIAPET_SET_DIR: $CHIAPET_SET_DIR"
echo "SAMPLE_DIRS: $SAMPLE_DIRS"

mkdir -p $CHIAPET_SET_DIR/log
LOG_FILE="log/preprocess_chiapet_set.log"

mkdir -p $CHIAPET_SET_DIR/peaks
echo "`date`: Run starting in $CHIAPET_SET_DIR" | tee $CHIAPET_SET_DIR/$LOG_FILE

rm -f $CHIAPET_SET_DIR/peaks/reads.bed
for SAMPLE_DIR in $SAMPLE_DIRS
do
    echo "`date`: Adding reads from $SAMPLE_DIR to reads.bed" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
    cat $SAMPLE_DIR/left.dedup.bed $SAMPLE_DIR/right.dedup.bed >> $CHIAPET_SET_DIR/peaks/reads.bed
done

echo "`date`: Calling peaks using reads.bed to define interaction anchor locations" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
which macs2
macs2 callpeak -t $CHIAPET_SET_DIR/peaks/reads.bed -f BED -n anchor --nomodel -p 0.01 --outdir $CHIAPET_SET_DIR/peaks
echo "`date`: Found `cat $CHIAPET_SET_DIR/peaks/anchor_peaks.narrowPeak | wc -l` peaks" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
bedtools merge -d $MERGE_GAP -i $CHIAPET_SET_DIR/peaks/anchor_peaks.narrowPeak > $CHIAPET_SET_DIR/peaks/anchor_peaks.merged.bed
echo "`date`: Merged peaks within ${MERGE_GAP}bp resulting in `cat $CHIAPET_SET_DIR/peaks/anchor_peaks.merged.bed | wc -l` anchors" | tee -a $CHIAPET_SET_DIR/$LOG_FILE


for SAMPLE_DIR in $SAMPLE_DIRS
do
    echo "`date`: Processing $SAMPLE_DIR" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
    cp $SAMPLE_DIR/bwa.log $CHIAPET_SET_DIR/log/`basename $SAMPLE_DIR`.bwa.log
    cp $SAMPLE_DIR/preprocess_fastq.log $CHIAPET_SET_DIR/log/`basename $SAMPLE_DIR`.preprocess_fastq.log    
    cp $SAMPLE_DIR/read_stats.txt $CHIAPET_SET_DIR/log/`basename $SAMPLE_DIR`.read_stats.txt
    echo "`date`:   Mapping `cat $SAMPLE_DIR/left.dedup.bed | wc -l` PETs to anchors by intersecting with peaks" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
    SAMPLE_ANCHOR_DIR=$CHIAPET_SET_DIR/peaks/`basename $SAMPLE_DIR`
    mkdir -p $SAMPLE_ANCHOR_DIR
    bedtools intersect -loj -a $SAMPLE_DIR/left.dedup.bed -b $CHIAPET_SET_DIR/peaks/anchor_peaks.merged.bed | awk '{print $5,$6,$7,$4}' OFS='\t' > $SAMPLE_ANCHOR_DIR/anchor1.bed
    bedtools intersect -loj -a $SAMPLE_DIR/right.dedup.bed -b $CHIAPET_SET_DIR/peaks/anchor_peaks.merged.bed | awk '{print $5,$6,$7,$4}' OFS='\t' > $SAMPLE_ANCHOR_DIR/anchor2.bed
    # Confirm reads match up between left and right
    cut -f 4 $SAMPLE_ANCHOR_DIR/anchor1.bed > $SAMPLE_ANCHOR_DIR/anchor1_names.txt
    cut -f 4 $SAMPLE_ANCHOR_DIR/anchor2.bed > $SAMPLE_ANCHOR_DIR/anchor2_names.txt
    DIFF=$(diff $SAMPLE_ANCHOR_DIR/anchor1_names.txt $SAMPLE_ANCHOR_DIR/anchor1_names.txt) 
    if [ "$DIFF" != "" ] 
    then
        echo "ERROR: Read names don't match between anchor1.bed and anchor2.bed" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
        exit 1
    fi
    # Create BEDPE
    paste $SAMPLE_ANCHOR_DIR/anchor1.bed $SAMPLE_ANCHOR_DIR/anchor2.bed | cut -f 1-3,5-8 > $SAMPLE_ANCHOR_DIR/anchor_interactions.tmp
    # Keep only interactions where both reads map to anchors
    awk '{if ($1 != "." && $4 != ".") print}' $SAMPLE_ANCHOR_DIR/anchor_interactions.tmp > $SAMPLE_ANCHOR_DIR/anchor_interactions.bedpe
    echo "`date`:   Wrote `cat $SAMPLE_ANCHOR_DIR/anchor_interactions.bedpe | wc -l` PETs where both reads map to anchors ($SAMPLE_ANCHOR_DIR/anchor_interactions.bedpe)" | tee -a $CHIAPET_SET_DIR/$LOG_FILE    
    cut -f1-6 $SAMPLE_ANCHOR_DIR/anchor_interactions.bedpe | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,".",$1}' > $SAMPLE_ANCHOR_DIR/loop_counts.bedpe
    cp $SAMPLE_ANCHOR_DIR/loop_counts.bedpe $CHIAPET_SET_DIR/`basename $SAMPLE_DIR`.loop_counts.bedpe
    echo "`date`:   Wrote loop PET counts to $CHIAPET_SET_DIR/`basename $SAMPLE_DIR`.loop_counts.bedpe" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
    SAME_ANCHOR_LOOPS=`cat $SAMPLE_ANCHOR_DIR/loop_counts.bedpe | awk '$1==$4 && $2==$5' | wc -l`
    DIFF_ANCHOR_LOOPS=`cat $SAMPLE_ANCHOR_DIR/loop_counts.bedpe | awk '$1!=$4 || $2!=$5' | wc -l`
    DIFF_ANCHOR_LOOPS_3PETS=`cat $SAMPLE_ANCHOR_DIR/loop_counts.bedpe | awk '($1!=$4 || $2!=$5) && $8>=3' | wc -l`
    DIFF_ANCHOR_INTRACHROMOSOMAL_LOOPS_3PETS=`cat $SAMPLE_ANCHOR_DIR/loop_counts.bedpe | awk '($1==$4) && ($2!=$5) && $8>=3' | wc -l`
    DIFF_ANCHOR_INTRACHROMOSOMAL_5kb_LOOPS_3PETS=`cat $SAMPLE_ANCHOR_DIR/loop_counts.bedpe | awk '($1==$4) && ($5-$3 >= 5000) && $8>=3' | wc -l`
    rm -fr $SAMPLE_ANCHOR_DIR
    echo "`date`:   Removed temporary directory $SAMPLE_ANCHOR_DIR" | tee -a $CHIAPET_SET_DIR/$LOG_FILE
    echo "`date`:   Invalid loops between the same anchor (i.e. left==right): $SAME_ANCHOR_LOOPS" | tee -a $CHIAPET_SET_DIR/$LOG_FILE    
    echo "`date`:   Valid loops between different anchors (i.e. left!=right): $DIFF_ANCHOR_LOOPS" | tee -a $CHIAPET_SET_DIR/$LOG_FILE        
    echo "`date`:   Valid loops between different anchors with 3+ supporting PETs: $DIFF_ANCHOR_LOOPS_3PETS" | tee -a $CHIAPET_SET_DIR/$LOG_FILE        
    echo "`date`:   Valid intrachromosomal loops between different anchors with 3+ supporting PETs: $DIFF_ANCHOR_INTRACHROMOSOMAL_LOOPS_3PETS" | tee -a $CHIAPET_SET_DIR/$LOG_FILE        
    echo "`date`:   Valid intrachromosomal 5kb+ loops with 3+ supporting PETs: $DIFF_ANCHOR_INTRACHROMOSOMAL_5kb_LOOPS_3PETS" | tee -a $CHIAPET_SET_DIR/$LOG_FILE        
done
