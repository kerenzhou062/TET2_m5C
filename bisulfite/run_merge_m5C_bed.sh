#!/bin/sh

function get_repel_bed {
  bed1=$1
  bed2=$2
  repel_bed=$3
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -v -s \
    > temp.bed

  awk 'BEGIN{FS="\t";OFS="\t";}
  ARGIND==1{
    if (FNR == 1) {
      header_line = $0;
    }
  }
  ARGIND==2{
    if (FNR == 1) {
      print header_line;
    }
    print $0;
  }' $bed1 temp.bed > ${repel_bed}

  rm -f temp.bed
}

function run_sbatch_merge_m5c_bed {
  sample_dir=$1
  input_dir=$2
  output_dir=$3
  type=$4
  if [[ "${type}" == "merge" ]]; then
    appendix=""
  else
    appendix=".${type}"
  fi
  ## make log directory
  if [[ ! -d "$OUTPUT_DIR" ]]; then
    mkdir -p $OUTPUT_DIR
  fi
  cd $OUTPUT_DIR
  ## loop sample files
  for sample_file in `find $sample_dir -type f -name "*.matrix" | grep -v "_multi"`;
  do
    base_name=$(basename $sample_file)
    prefix=${base_name%%.matrix}
    ## remove intron
    exon_bed=${INPUT_DIR}/${prefix}.exon${appendix}.bed
    intron_bed=${INPUT_DIR}/${prefix}.intron${appendix}.bed
    multi_exon_bed=${INPUT_DIR}/${prefix}_multi.exon${appendix}.bed
    multi_intron_bed=${INPUT_DIR}/${prefix}_multi.intron${appendix}.bed
    ## unique in intron
    bedtools intersect -a $intron_bed -b $exon_bed -s -v > intron.bed
    ## unique in multi exon
    cat $exon_bed intron.bed | bedtools intersect -a $multi_exon_bed -b stdin -s  -v > multi_exon.bed
    ## unique in multi intron
    cat $exon_bed intron.bed multi_exon.bed | bedtools intersect -a $multi_intron_bed -b stdin -s -v > multi_intron.bed
    cat $exon_bed intron.bed multi_exon.bed multi_intron.bed > ${prefix}.merge.bed.tmp
    ## sort bed files
    awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR==1){print $0}}' ${prefix}.merge.bed.tmp > header.tmp
    awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR>1){print $0}}' ${prefix}.merge.bed.tmp | sort -t$'\t' -k1,1 -k2,2 > bed.tmp
    cat header.tmp bed.tmp | \
      awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR==1){print $0}else{$4=$4"="FNR; print $0}}' > ${prefix}.merge${appendix}.bed
    ## delete temp files
    rm -f intron.bed multi_exon.bed multi_intron.bed *.tmp
  done
}

function get_hypo_bed {
  bed1=$1
  bed2=$2
  output=$3
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -v -s | \
    awk 'BEGIN{FS="\t";OFS="\t";}{print $0}' > bed1.uniq.tmp
  bedtools intersect \
    -a $bed1 \
    -b $bed2 -wa -wb -s | \
    awk 'BEGIN{FS="\t";OFS="\t";cutoff=0.05;}{
      name=labela"_"labelb"_"FNR;
      diff = $5-$16;
      if (diff > cutoff) {
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11;
      }
    }' > intersect.tmp
  head -n 1 $bed1 > header.tmp
  #with random values
  cat header.tmp bed1.uniq.tmp intersect.tmp > $output
  rm -f *.tmp
}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

## for caRNA-batch1
INPUT_DIR="$BASE/m5C_site/caRNA/batch1/bed"
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed"
LOG_DIR="$BASE/log/caRNA/call_m5c/batch1"

## 
SAMPLE_DIR="$BASE/sample/caRNA/RNA-m5C/batch1/exon/"
run_sbatch_merge_m5c_bed $SAMPLE_DIR $INPUT_DIR $OUTPUT_DIR 'merge'

cd $OUTPUT_DIR
get_hypo_bed PDX148_EV_bisulfite.merge.bed PDX148_TET2-bisulfite.merge.bed PDX148_m5C_site_repress_by_TET2.bed
