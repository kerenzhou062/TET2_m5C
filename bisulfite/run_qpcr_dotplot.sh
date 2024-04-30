#!/bin/sh

function creat_folder {
  if [ ! -d "$1" ]; then
    mkdir -p "$1"
  fi
}


BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

OUTPUT_DIR="$BASE/analysis/qpcr_dotplot"

creat_folder $OUTPUT_DIR
cd $OUTPUT_DIR

$SCRIPTS/get_transform_qpcr.py ./ THP1_KOOE.combine.txt
$SCRIPTS/draw_dotplot.r
