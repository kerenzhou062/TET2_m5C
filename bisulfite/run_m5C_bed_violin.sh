#!/bin/sh

function creat_folder {
  if [ ! -d "$1" ]; then
    mkdir -p "$1"
  fi
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
  echo -e "id\tx\ty\tstatus" > header.txt
  cat header.txt bed1.uniq.txt bed2.uniq.txt intersect.txt > $output
  rm -f *.txt
}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

## for caRNA-batch1
INPUT_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed"
OUTPUT_DIR="$BASE/analysis/violin/batch1"

creat_folder $OUTPUT_DIR
cd $OUTPUT_DIR

get_scatter_matrix $INPUT_DIR/PDX148_EV_bisulfite.merge.bed $INPUT_DIR/PDX148_TET2-bisulfite.merge.bed PDX148_EV PDX148_TET2 PDX148_EV_TET2.scatter.matrix

creat_folder $OUTPUT_DIR
cd $OUTPUT_DIR

ANNO_DIR="$BASE/analysis/scatter/batch1"
## label scatter with label
awk 'BEGIN{FS="\t";OFS="\t";}
ARGIND==1{
  if (FNR > 1) {
    if ($21 == "ERV1" || $21 == "ERVL" || $21 == "ERVK" || $21 == "ERVL-MaLR") {
      type = "LTR";
    }else{
      type = $21;
    }
    arr[$4] = type;
  }
}
ARGIND==2{
  if (FNR >1) {
    type = arr[$1];
    if (type == "LTR" || type == "L1" || type == "L2") {
      print type, "EV",  $2;
      print type, "TET2", $3;
    }
  }else{
    print "type", "condition", "level";
  }
}' $ANNO_DIR/PDX148_EV_TET2.scatter.anno.txt PDX148_EV_TET2.scatter.matrix >  PDX148_EV_TET2.violin.matrix

$SCRIPTS/draw_violin_plot.r
