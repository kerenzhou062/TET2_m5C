#!/bin/bash

function Run_meme_motif {
  ## delete and remake
  INPUT_DIR=$1
  RESULT_DIR=$2
  if [ ! -d "$RESULT_DIR" ]; then
   mkdir -p $RESULT_DIR
  fi
  echo "" > $RESULT_DIR/command.log
  ## annotation
  for bed in $BED_DIR/*.bed; do
    INPUT_BED=$bed
    baseName=${INPUT_BED%%.bed}
    prefix=$(basename $baseName)
    BED_PATH=`realpath $INPUT_BED`
    ## run
    echo "bsub -P ENCORI -q rhel8_standard -M 20000 -n 10 \
      -eo $LOG_DIR/$prefix.err \
      -oo $LOG_DIR/$prefix.log \
      -J homer_motif_${prefix} \
      \"bash $SCRIPTS/sbatch_m5C_homer_motif.sh \
      $BED_PATH \
      $RESULT_DIR \
      $prefix\"" >> $RESULT_DIR/command.log
    ## submit jobs
    bsub -P ENCORI -q rhel8_standard -M 20000 -n 10 \
      -eo $LOG_DIR/$prefix.meme.err \
      -oo $LOG_DIR/$prefix.meme.log \
      -J meme_motif_${prefix} \
      "bash $SCRIPTS/sbatch_m5C_meme_motif.sh \
      $GENOME_FA \
      $GENOME_SIZE \
      $BED_PATH \
      $RESULT_DIR \
      $prefix"
  done
}

### running analysis ###
root=/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou
GENOME_FA=$root/public/genome/fasta/hg38/hg38.fa
GENOME_SIZE=$root/public/genome/size/hg38.chrom.sizes

## annotation
BASE=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
SCRIPTS=$BASE/scripts
BED_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed"
ANALYSIS_DIR="$BASE/analysis/motif_meme/caRNA/batch1"
LOG_DIR="$BASE/log/motif/caRNA/batch1"

rm -rf $LOG_DIR
mkdir -p $LOG_DIR

Run_meme_motif $BED_DIR $ANALYSIS_DIR
