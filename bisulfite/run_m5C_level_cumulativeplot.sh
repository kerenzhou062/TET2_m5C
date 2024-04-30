#!/bin/bash

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data`
LOG_DIR="$BASE/log"
SCRIPTS="$BASE/scripts"

BED_DIR="$BASE/m5C_site/caRNA/batch1/merge_bed/"
OUTPUT_DIR="$BASE/analysis/comulative_plot"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

## get pure TE
awk 'BEGIN{FS="\t";OFS="\t";} {if (FNR > 1) print $5}' $BED_DIR/PDX148_EV_bisulfite.merge.bed > PDX148_EV_bisulfite.m5c_level.txt
awk 'BEGIN{FS="\t";OFS="\t";} {if (FNR > 1) print $5}' $BED_DIR/PDX148_TET2-bisulfite.merge.bed > PDX148_TET2_bisulfite.m5c_level.txt

$SCRIPTS/cumulativePlot.pl -header false -interval 0.01 \
  --input PDX148_EV_bisulfite.m5c_level.txt \
    PDX148_TET2_bisulfite.m5c_level.txt \
  -o PDX148_EV_TET2_bisulfite.cumulativePlot.txt

sed -i '1i FC\tPDX148_EV\tPDX148_TET2' PDX148_EV_TET2_bisulfite.cumulativePlot.txt

$SCRIPTS/draw_cumulative_plot.r \
  $OUTPUT_DIR \
  PDX148_EV_TET2_bisulfite.cumulativePlot.txt \
  PDX148_EV_TET2_bisulfite.cumulativePlot.pdf "TRUE" \
  PDX148_EV_bisulfite.m5c_level.txt \
  PDX148_TET2_bisulfite.m5c_level.txt \
  PDX148_TET2
