#!/bin/sh

function run_sbatch_cutadapt {
  matrix_file=$1
  input_fastq_folder=$2
  output_fastq_folder=$3
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
    bsub -P cutadapt -q rhel8_standard -n ${THREAD} -M 50000 -W "240:00" \
      -e $log_folder/$exp_prefix.cutadapt.err -o $log_folder/$exp_prefix.cutadapt.log \
      -J $exp_prefix.cutadapt \
      "$SCRIPTS/sbatch_cutadapt.sh $input_fastq_folder $output_fastq_folder $exp_prefix $THREAD"
  done < $matrix_file

}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=5


## for caRNA-batch1
MATRIX_DIR="$BASE/sample/caRNA"
LOG_DIR="$BASE/log/caRNA/cutadapt/batch1"
FASTQ_DIR="$BASE/fastq/caRNA/batch1"
run_sbatch_cutadapt $MATRIX_DIR/sample.batch1.matrix $FASTQ_DIR/original $FASTQ_DIR $LOG_DIR

