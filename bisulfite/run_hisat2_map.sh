#!/bin/sh

function run_sbatch_hisat2_map {
  matrix_file=$1
  input_fastq_folder=$2
  align_folder=$3
  log_folder=$4
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
      bsub -P hisat2_map -q rhel8_standard -n 16 -M 150000 -W "240:00" \
        -eo $log_folder/$exp_prefix.hisat2_map.err -oo $log_folder/$exp_prefix.hisat2_map.log \
        -J $exp_prefix.hisat2_map \
        "$SCRIPTS/sbatch_hisat2_map.sh $input_fastq_folder $align_folder/$exp_prefix $exp_prefix"
    fi
  done < $matrix_file

}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=5

## for caRNA-batch1
MATRIX_DIR="$BASE/sample/caRNA"
LOG_DIR="$BASE/log/caRNA/alignment/batch1"
FASTQ_DIR="$BASE/fastq/caRNA/batch1"
ALIGN_DIR="$BASE/alignment/caRNA/batch1"

run_sbatch_hisat2_map $MATRIX_DIR/sample.batch1.matrix $FASTQ_DIR $ALIGN_DIR $LOG_DIR

## for caRNA-batch2
MATRIX_DIR="$BASE/sample/caRNA"
LOG_DIR="$BASE/log/caRNA/alignment/batch2"
FASTQ_DIR="$BASE/fastq/caRNA/batch2"
ALIGN_DIR="$BASE/alignment/caRNA/batch2"

run_sbatch_hisat2_map $MATRIX_DIR/sample.batch2.matrix $FASTQ_DIR $ALIGN_DIR $LOG_DIR
