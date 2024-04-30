#!/bin/sh

function run_sbatch_pipeup_intron {
  matrix_file=$1
  align_folder=$2
  log_folder=$3
  ## make log directory
  if [[ ! -d "$log_folder" ]]; then
    mkdir -p $log_folder
  fi
  ## loop file by line
  i=1
  while IFS=$'\t' read -r -a matrixArr
  do
    test $i -eq 1 && ((i=i+1)) && continue
    exp_prefix="${matrixArr[0]}"
    tech="${matrixArr[7]}"
    if [ "$tech" == "bisulfite" ]; then
      bsub -P pipeup_intron -q rhel8_standard -n 11 -M 80000 -W "240:00" \
        -eo $log_folder/$exp_prefix.pipeup_intron.err -oo $log_folder/$exp_prefix.pipeup_intron.log \
        -J $exp_prefix.pipeup_intron \
        "$SCRIPTS/sbatch_pipeup_intron.sh $align_folder/$exp_prefix $exp_prefix"
    fi
  done < $matrix_file

}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=5

## for caRNA-batch1
MATRIX_DIR="$BASE/sample/caRNA"
LOG_DIR="$BASE/log/caRNA/alignment/batch1"
ALIGN_DIR="$BASE/alignment/caRNA/batch1"

run_sbatch_pipeup_intron $MATRIX_DIR/sample.batch1.matrix $ALIGN_DIR $LOG_DIR
