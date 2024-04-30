#!/bin/sh
BASE="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data"

OUTPUT_DIR="$BASE/analysis/DESeq2/batch1"
# GENOME
FULL_BED="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/annotation/hg38.gencode.v37.annotation.caRNA.anno.bed12"

if [[ ! -d "$OUTPUT_DIR" ]]; then
  mkdir -p $OUTPUT_DIR
fi

cd $OUTPUT_DIR

##get counts of genes
buildExpMatrix.py -input $BASE/alignment_unconvert/caRNA/batch1 -identity genes -grepKept 'PDX148' -grepExpel 'IFN' -abundance counts \
  -output genes.counts.matrix.tmp

##get TPM of isoforms
buildExpMatrix.py -input $BASE/alignment_unconvert/caRNA/batch1 -identity isoforms -grepKept 'PDX148' -grepExpel 'IFN' -abundance counts \
  -output isoforms.counts.matrix.tmp
  
awk 'BEGIN{FS="\t";OFS="\t";}{if (FNR==1){gsub(/-/, "_", $0)} print $0}' genes.counts.matrix.tmp > genes.counts.matrix
awk 'BEGIN{FS="\t";OFS="\t";}{if (FNR==1){gsub(/-/, "_", $0)} print $0}' isoforms.counts.matrix.tmp > isoforms.counts.matrix
rm -f *.tmp
## run DESeq2 on samples
## differential gene expression analysis
function RunDESeq2 {
  matrix=$1
  counts=$2
  control=$3
  treat=$4
  design=$5
  prefix=$6
  directory=$7
  batch=$8
  args=""
  output_prefix="$prefix.unbatch"
  ## if run with batch effect correction
  if [[ "$batch" == "true" ]]; then
    args=" --batchMethod spikeins "
    output_prefix="$prefix.batch"
  fi

  DESeq2Gene.R --sampleMtx $matrix --counts $counts \
    --control $control --treat $treat --design $design -f "0" $args --relevel "$design:$control" \
    --subset --test 'Wald' --pval 0.05 --adjp 1 --shrink 'apeglm' --prefix $output_prefix --output $directory
}

#
COUNT_MATRIX="$OUTPUT_DIR/genes.counts.matrix"
SAMPLE_MATRIX="$BASE/sample/caRNA/sample.batch1.DESeq2.matrix"

RunDESeq2 $SAMPLE_MATRIX $COUNT_MATRIX "Ctrl_EV" "Ctrl_TET2" "condition" "PDX148_TET2_vs_EV_gene" ./ False
