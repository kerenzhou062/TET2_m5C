#!/bin/sh

function creat_folder {
  if [ ! -d "$1" ]; then
    mkdir -p "$1"
  fi
}

function renameGeneType {
  input=$1
  $SCRIPTS/refine_m5c_target_geneType.py \
  --protein $PROTEIN_CLASS \
  --anno $CUSTOM_GENE_TYPE \
  --rnasub $RNA_SUBLOCATION \
  --input $input \
  --output $input.tmp
  mv $input.tmp $input
}

function get_scatter_matrix {
  bed1=$1
  bed2=$2
  label1=$3
  label2=$4
  output=$5
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -v -s | \
    awk -v label="$label1" 'BEGIN{FS="\t";OFS="\t";}{name=label"_"FNR;print name, $5, 0, "Hypo"}' > bed1.uniq.txt
  bedtools intersect \
    -a $bed2 \
    -b $bed1 -v -s | \
    awk -v label="$label2" 'BEGIN{FS="\t";OFS="\t";}{name=label"_"FNR;print name, 0, $5, "Hyper"}' > bed2.uniq.txt
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -wa -wb -s | \
    awk -v labela="$label1" -v labelb="$label2" 'BEGIN{FS="\t";OFS="\t";cutoff=0.05;}{
      name=labela"_"labelb"_"FNR;
      diff = $5-$16;
      if (diff > cutoff) {
        label="Hypo";
      }else if (diff < -cutoff) {
        label="Hyper";
      }else{
        label="nochange";
      }
      print name, $5, $16, label;
    }' > intersect.txt
  #with random values
  $SCRIPTS/get_random_scatter_value.py bed1.uniq.txt 2 0.03 100 > bed1.uniq.random.txt
  $SCRIPTS/get_random_scatter_value.py bed2.uniq.txt 1 0.03 200 > bed2.uniq.random.txt
  echo -e "id\tx\ty\tstatus" > header.txt
  cat header.txt bed1.uniq.random.txt bed2.uniq.random.txt intersect.txt > $output
  rm -f *.txt
}

function get_scatter_bed {
  bed1=$1
  bed2=$2
  label1=$3
  label2=$4
  output=$5
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -v -s | \
    awk -v label="$label1" 'BEGIN{FS="\t";OFS="\t";}{name=label"_"FNR;$4=name; print $0}' > bed1.uniq.bed
  bedtools intersect \
    -a $bed2 \
    -b $bed1 -v -s | \
    awk -v label="$label2" 'BEGIN{FS="\t";OFS="\t";}{name=label"_"FNR;$4=name; print $0}' > bed2.uniq.bed
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -wa -wb -s | \
    awk -v labela="$label1" -v labelb="$label2" 'BEGIN{FS="\t";OFS="\t";cutoff=0.05;}{
      name=labela"_"labelb"_"FNR;
      diff = $5-$16;
      if (diff > cutoff) {
        label="Hypo";
      }else if (diff < -cutoff) {
        label="Hyper";
      }else{
        label="nochange";
      }
      $4=name;
      print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11;
    }' > intersect.bed
  #with random values
  head -n 1 $bed1 > header.bed
  cat header.bed bed1.uniq.bed bed2.uniq.bed intersect.bed > $output
  rm -f header.bed bed1.uniq.bed bed2.uniq.bed intersect.bed
}

function get_scatter_label {
  scatter_matrix=$1
  annores=$2
  label=$3
  output=$4
  ## label="L1M1|L1HS|LTR2B|LTR2C|LTR12C|LTR5|LTR13A"
  awk -v label="$label" 'BEGIN{FS="\t";OFS="\t";split(label, label_arr, "|");for(i in label_arr){each=label_arr[i];label_arr[each]=each;}}
  ARGIND==1{
    if ($14 in label_arr) {
      anno_arr[$4] = $14;
    }
  }
  ARGIND==2{
    if (FNR == 1) {
      print $0, "gene";
    }else{
      if ($1 in anno_arr) {
        if ($4 == "Hypo") {
          print $0, anno_arr[$1];
        }else{
          print $0, "NA";
        }
      }else{
        print $0, "NA";
      }
    }
  }' $annores $scatter_matrix > $output
}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

## for caRNA-batch1
INPUT_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed"
OUTPUT_DIR="$BASE/analysis/scatter/batch1"

creat_folder $OUTPUT_DIR
cd $OUTPUT_DIR

get_scatter_matrix $INPUT_DIR/PDX148_EV_bisulfite.merge.bed $INPUT_DIR/PDX148_TET2-bisulfite.merge.bed PDX148_EV PDX148_TET2 PDX148_EV_TET2.scatter.matrix
get_scatter_bed $INPUT_DIR/PDX148_EV_bisulfite.merge.bed $INPUT_DIR/PDX148_TET2-bisulfite.merge.bed PDX148_EV PDX148_TET2 PDX148_EV_TET2.scatter.bed

## genome
ROOT=/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou
PUBLIC_BASE="$ROOT/public/genome"
GENOME_SIZE="$PUBLIC_BASE/size/hg38.chrom.sizes"
NCBI_GENE_INFO="$PUBLIC_BASE/annotation/hg38/hg38.ncbi.gene_info"

## transcript expression
EXP_MTX=$ROOT/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/analysis/DESeq2/caRNA_tx/caRNA_tx.TPM.txt

## annotation
BASE=$ROOT/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
MAIN_GENE_TYPE="$BASE/public/geneType.hg38.v37.txt"
CUSTOM_GENE_TYPE="$BASE/public/gene_type.txt"
PROTEIN_CLASS="$BASE/public/protein_class.proteinatlas.txt"
RNA_SUBLOCATION="$BASE/public/All_RNA_subcellular_localization_data.RNALocate.txt"

SCRIPTS=$BASE/scripts
CAR_RNA_BED="$ROOT/project/chen_yang_lab/ying_qing/tet1/m5c_seq/inhouse_data/annotation/hg38.gencode.v37.annotation.caRNA.anno.bed12"
peakAnnotate.py -input PDX148_EV_TET2.scatter.bed \
  --keepAll -codon 25,25 \
  -by 'gene' \
  -priF "start_codon,stop_codon,5' UTR,CDS,3' UTR,exon,intron" \
  -priG "protein_coding,IG_gene,TR_gene,lncRNA,pseudogene,miRNA,snoRNA,rRNA,tRNA,snRNA,scRNA,other_sRNA,paRNA,eRNA" \
  --keepName -peakType "bed6" -mode RNA -gsize $GENOME_SIZE \
  -anno $CAR_RNA_BED -geneClassFile $MAIN_GENE_TYPE \
  -ncbiGeneInfo $NCBI_GENE_INFO \
  -output PDX148_EV_TET2.scatter.anno.txt

renameGeneType PDX148_EV_TET2.scatter.anno.txt

## label scatter with label
label="L1HS|L1M1|LTR5|LTR13A"
get_scatter_label PDX148_EV_TET2.scatter.matrix PDX148_EV_TET2.scatter.anno.txt  $label PDX148_EV_TET2.scatter.label.matrix
#$SCRIPTS/draw_scatter_plot.r
