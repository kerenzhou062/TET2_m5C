#!/bin/bash

function RunPeakGeneLink {
  ## delete and remake
  INPUT_DIR=$1
  RESULT_DIR=$2
  STATS_DIR=$3
  EXP_MTX=$4
  rm -rf $ANNODATION_DIR
  mkdir -p $ANNODATION_DIR
  rm -rf $STATS_DIR
  mkdir -p $STATS_DIR
  #
  cd $RESULT_DIR
  #
  find "$INPUT_DIR" -type f -name "*.bed" | xargs -I {} cp {} ./

  echo "" > $ANALYSIS_DIR/command.log
  ## annotation
  for bed in *.bed; do
    INPUT_BED=$bed
    baseName=${INPUT_BED%%.bed}
    prefix=$(basename $baseName)
    BED_PATH=`realpath $INPUT_BED`
    #command="AnnoBed $INPUT_BED $RESULT_DIR $STATS_DIR"
    echo "bsub -P ENCORI -q rhel8_standard -M 10000 -n 1 \
      -eo $LOG_DIR/$prefix.err \
      -oo $LOG_DIR/$prefix.log \
      -J RBP_target_${prefix} \
      \"bash $SCRIPTS/sbatch_m5C_annotation.sh \
      $BED_PATH \
      $RESULT_DIR \
      $STATS_DIR \
      $MAIN_GENE_TYPE \
      $NCBI_GENE_INFO \
      $GENOME_SIZE \
      $CAR_RNA_BED \
      $RNA_SUBLOCATION \
      $PROTEIN_CLASS \
      $CUSTOM_GENE_TYPE \
      $SCRIPTS \
      $EXP_MTX \"" >> $ANALYSIS_DIR/command.log
    ## run
    bsub -P ENCORI -q rhel8_standard -M 10000 -n 1 \
      -eo $LOG_DIR/$prefix.err \
      -oo $LOG_DIR/$prefix.log \
      -J RBP_target_${prefix} \
      "bash $SCRIPTS/sbatch_m5C_annotation.sh \
      $BED_PATH \
      $RESULT_DIR \
      $STATS_DIR \
      $MAIN_GENE_TYPE \
      $NCBI_GENE_INFO \
      $GENOME_SIZE \
      $CAR_RNA_BED \
      $RNA_SUBLOCATION \
      $PROTEIN_CLASS \
      $CUSTOM_GENE_TYPE \
      $SCRIPTS \
      $EXP_MTX"
  done
}

### running analysis ###
root=/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou

## genome
PUBLIC_BASE="$root/public/genome"
GENOME_SIZE="$PUBLIC_BASE/size/hg38.chrom.sizes"
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"

## transcript expression
EXP_MTX=$root/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/analysis/DESeq2/caRNA_tx/caRNA_tx.TPM.txt

## annotation
BASE=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
MAIN_GENE_TYPE="$BASE/public/geneType.hg38.v37.txt"
CUSTOM_GENE_TYPE="$BASE/public/gene_type.txt"
PROTEIN_CLASS="$BASE/public/protein_class.proteinatlas.txt"
RNA_SUBLOCATION="$BASE/public/All_RNA_subcellular_localization_data.RNALocate.txt"

SCRIPTS=$BASE/scripts
CAR_RNA_BED="$root/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/annotation/hg38.gencode.v37.annotation.caRNA.anno.bed12"
LOG_DIR="$BASE/log/annotation"

## for batch1
BED_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed"
ANALYSIS_DIR="$BASE/m5C_annotation/caRNA/batch1"
ANNODATION_DIR="$ANALYSIS_DIR/annotation"
STATS_DIR="$ANALYSIS_DIR/stats"

rm -rf $LOG_DIR
mkdir -p $LOG_DIR

RunPeakGeneLink "$BED_DIR" $ANNODATION_DIR $STATS_DIR $EXP_MTX

## for batch2
BED_DIR="$BASE/m5C_site/caRNA/batch2/merge_bed"
ANALYSIS_DIR="$BASE/m5C_annotation/caRNA/batch2"
ANNODATION_DIR="$ANALYSIS_DIR/annotation"
STATS_DIR="$ANALYSIS_DIR/stats"

rm -rf $LOG_DIR
mkdir -p $LOG_DIR

RunPeakGeneLink "$BED_DIR" $ANNODATION_DIR $STATS_DIR $EXP_MTX
