#!/bin/sh
#BSUB -n 39
#BSUB -R "rusage[mem=5GB]"
#BSUB -W 48:00
#BSUB -J run_get_caRNA_expression.gene
#BSUB -P run_get_caRNA_expression.gene
#BSUB -N
#BSUB -oo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/log/run_get_caRNA_expression.gene.log
#BSUB -eo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/log/run_get_caRNA_expression.gene.err

function link_bam {
  input_dir=$1
  find $input_dir -type f -name "*.genome.sorted.bam" | xargs -I {} ln -s {} ./
  find $input_dir -type f -name "*.genome.sorted.bam.bai" | xargs -I {} ln -s {} ./
}

function unlink_bam {
  for i in *.bam;
  do
    if [[ -L ${i} ]]; then
      unlink $i
    fi
  done
  for i in *.bam.bai;
  do
    if [[ -L ${i} ]]; then
      unlink $i
    fi
  done
}

function run_featureCounts {
  gtf=$1
  att=$2
  output_dir=$3
  output=$4
  if [[ $att == "gene_id" ]]; then
    feature="gene"
  else
    feature="transcript"
  fi
  featureCounts -a $gtf -F GTF \
    -g $att -t $feature \
    -s 2 -T 38 -p \
    --fracOverlap 0 -C -M -O \
    --fraction --countReadPairs \
    --donotsort \
    -o $output \
    $output_dir/*.bam
}

function run_get_exp {
  input_dir=$1
  output_dir=$2
  prefix=$3
  gtf=$4
  att=$5
  unlink_bam $output_dir
  link_bam $input_dir
  run_featureCounts $gtf $att $output_dir $prefix.raw.txt
  ## normalize gene counts
  $SCRIPTS_DIR/format_featureCounts.py \
    --input $output_dir/$prefix.raw.txt \
    --method counts \
    --output $output_dir/$prefix.counts.txt
  $SCRIPTS_DIR/format_featureCounts.py \
    --input $output_dir/$prefix.raw.txt \
    --method TPM \
    --output $output_dir/$prefix.TPM.txt
  $SCRIPTS_DIR/format_featureCounts.py \
    --input $output_dir/$prefix.raw.txt \
    --method FPKM \
    --output $output_dir/$prefix.FPKM.txt
  ## unlink soft links
  unlink_bam $output_dir
}

## run testing
BASE="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data"
SCRIPTS_DIR="$BASE/scripts"
# GENOME
PUBLIC_BASE="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/public/genome"
GENETYPE_FILE="$PUBLIC_BASE/annotation/hg38/geneType.hg38.txt"
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"
FULL_BED="$PUBLIC_BASE/annotation/hg38/v37/gencode.v37.annotation.anno.bed12"

# evaluate sample distances

## caRNA
BAM_DIR="$BASE/alignment/caRNA/batch2"
OUTPUT_DIR="$BASE/analysis/DESeq2/caRNA_gene"
GTF="$BASE/annotation/hg38.gencode.v37.annotation.caRNA.gtf"

rm -rf $OUTPUT_DIR
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

run_get_exp $BAM_DIR $OUTPUT_DIR caRNA_gene $GTF "gene_id"
