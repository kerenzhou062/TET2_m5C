#!/bin/bash

function append_exp {
  anno=$1
  exp_mtx=$2
  echo $anno
  echo $exp_mtx
  awk 'BEGIN{FS="\t";OFS="\t";}
  ARGIND==1{
    if (FNR >1) {
      tpm = ($2 + $3) / 2;
      arr[$1] = tpm;
    }
  }
  ARGIND==2{
    if (FNR == 1) {
      print $0, "transcript_TPM";
    }else{
      print $0, arr[$19];
    }
  }' $exp_mtx $anno > $anno.tmp
  mv $anno.tmp $anno
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

function getGeneTypeExp {
  input=$1
  prefix=${2}.geneTypeExp
  $SCRIPTS/determin_m5c_target_geneTypeExp.py \
    --input $input \
    --output $prefix.txt
}

function plotRnaSublocation {
  input=$1
  prefix=${2}.rnaSublocation
  $SCRIPTS/determin_m5c_target_sublocation.py \
    --input $input \
    --output $prefix.txt
  drawPieChart.py --input $prefix.txt \
  --label "type" --value "count" --keep --opacity 0.9 \
  --order "Nucleolus,Nucleoplasm,Nucleus,Nuclear,Chromatin,Cytoplasm,Insoluble cytoplasm,Cytosol,Ribosome-free cytosol,Ribosome,ER-bound polysome,Cytosolic polysome,Endoplasmic reticulum,Mitochondrion,Microvesicle,Membrane,Exosome,Extracellular vesicle,Other,Unkown" \
  --color "#499894,#86BCB6,#F28E2B,#FFBE7D,#59A14F,#8CD17D,#B6992D,#F1CE63,#9D7660,#D7B5A6,#E15759,#FF9D9A,#A0CBE8,#4E79A7,#D37295,#FABFD2,#B07AA1,#D4A6C8,#79706E,#BAB0AC" \
  --title "Rna sublocation statistics" \
  --output $prefix.pdf > /dev/null 2>&1
  #pdftoppm $prefix.pdf $prefix \
  #  -rx 600 -ry 600 -png > /dev/null 2>&1
}

function plotGeneType {
  input=$1
  prefix=${2}.geneType
  $SCRIPTS/determin_m5c_target_geneType.py \
    --input $input \
    --output $prefix.txt
  drawPieChart.py --input $prefix.txt \
  --label "type" --value "count" --keep --opacity 0.8 \
  --order "CAC,CD,ENZYME,GPCRs,Kinases,Metabolic,NRs,RAS,Ribosomal,TF,Transporters,Peptidases,Plasma,pMembrane,pSecreted,Other_protein,miRNA,lncRNA,pseudogene,rRNA,snoRNA,snRNA,scRNA,tRNA,other_sRNA,intergenic" \
  --color "#ffa500,#b54871,#4e8c8c,#9acd32,#44a944,#de946d,#6a6dc1,#e56de5,#8e5ff0,#bc8f8f,#4682b4,#baba27,#8fbc8f,#e63636,#e4813a,#616161,#943dd3,#2dc477,#f0e68c,#4e64ee,#1e90ff,#f4dd1b,#dda0dd,#e2627c,#f366b2,#999797" \
  --title "Gene type statistics" \
  --output $prefix.pdf > /dev/null 2>&1
  #pdftoppm $prefix.pdf $prefix \
  #  -rx 600 -ry 600 -png > /dev/null 2>&1
}

function plotTxType {
  input=$1
  prefix=${2}.tx_type
  $SCRIPTS/determin_m5c_target_tx_type.py \
    --input $input \
    --output $prefix.all.txt
  $SCRIPTS/determin_m5c_target_tx_type.repeat.py \
    --input $input \
    --output $prefix.repeat.txt
  drawPieChart.py --input $prefix.repeat.txt \
  --label "type" --value "count" --keep --opacity 0.8 \
  --order "Alu,CR1,ERV1,ERVL,ERVK,ERVL-MaLR,MIR,hAT-Charlie,hAT-Tip100,hAT-Blackjack,L1,L2,RTE-X,TcMar-Tigger,Other repeat,snRNA,tRNA" \
  --color "#E15759,#FF9D9A,#9D7660,#D7B5A6,#D37295,#FABFD2,#B07AA1,#D4A6C8,#4E79A7,#A0CBE8,#F28E2B,#FFBE7D,#59A14F,#8CD17D,#B6992D,#F1CE63,#499894" \
  --title "Tx type statistics" \
  --output $prefix.repeat.pdf > /dev/null 2>&1
  #pdftoppm $prefix.pdf $prefix \
  #  -rx 600 -ry 600 -png > /dev/null 2>&1
}

function plotRepeatInRegion {
  input=$1
  prefix=${2}.repeat
  $SCRIPTS/determin_m5c_target_repeat_in_region.py \
    --input $input \
    --output $prefix.txt
  drawSunburst.py --input $prefix.txt \
  --label "type,repeat" --value "count" --keep --opacity 0.8 \
  --order "Exon,5-UTR,start_codon,CDS,stop_codon,3-UTR,Coding exon,Non-coding exon,Intron,Alu,ERV1,ERVL-MaLR,ERVL,hAT-Charlie,CR1,L1,L2,MIR,TcMar-Tigger,Other repeat,intergenic" \
  --color "#58595B,#4E79A7,#A0CBE8,#F28E2B,#FFBE7D,#59A14F,#8CD17D,#B6992D,#F1CE63,#499894,#86BCB6,#E15759,#FF9D9A,#79706E,#BAB0AC,#D37295,#FABFD2,#B07AA1,#D4A6C8,#9D7660,#D7B5A6" \
  --title "Repeat in region statistics" \
  --output $prefix.pdf > /dev/null 2>&1
  #pdftoppm $prefix.pdf $prefix \
  #  -rx 600 -ry 600 -png > /dev/null 2>&1
}

########### annotation ##########
function AnnoBed {
  inputBed=$1
  output_dir=$2
  stats_dir=$3
  exp_mtx=$4
  baseName=${inputBed%%.bed}
  rbp_label=$(basename $baseName)
  ## assign only one possible target gene to a peak
  # annotate peaks
  output=$output_dir/${rbp_label}.peakAll.anno.txt
  peakAnnotate.py -input $inputBed \
    --keepAll -codon 25,25 \
    -by 'gene' \
    -priF "start_codon,stop_codon,5' UTR,CDS,3' UTR,exon,intron" \
    -priG "protein_coding,IG_gene,TR_gene,lncRNA,pseudogene,miRNA,snoRNA,rRNA,tRNA,snRNA,scRNA,other_sRNA,paRNA,eRNA" \
    --keepName -peakType "bed6" -mode RNA -gsize $GENOME_SIZE \
    -anno $CAR_RNA_BED -geneClassFile $MAIN_GENE_TYPE \
    -ncbiGeneInfo $NCBI_GENE_INFO \
    -output $output
  renameGeneType $output
  #plotGeneType $output $stats_dir/${rbp_label}.peakAll
  plotTxType $output $stats_dir/${rbp_label}.peakAll
  #plotRepeatInRegion $output $stats_dir/${rbp_label}.peakAll
  #plotRnaSublocation $output $stats_dir/${rbp_label}.peakAll
  append_exp $output $exp_mtx
}

inputBed=$1
output_dir=$2
stats_dir=$3
MAIN_GENE_TYPE=$4
NCBI_GENE_INFO=$5
GENOME_SIZE=$6
CAR_RNA_BED=$7
RNA_SUBLOCATION=${8}
PROTEIN_CLASS=${9}
CUSTOM_GENE_TYPE=${10}
SCRIPTS=${11}
exp_mtx=${12}

AnnoBed $inputBed $output_dir $stats_dir $exp_mtx
