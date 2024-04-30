#!/bin/sh

function create_folder {
  if [ ! -d "$1" ]; then
      mkdir -p "$1"
  fi
}

function run_sbatch_call_m5c {
  sample_dir=$1
  output_dir=$2
  log_folder=$3
  type=$4
  ## make log directory
  if [[ ! -d "$log_folder" ]]; then
    mkdir -p $log_folder
  fi
  ## loop sample files
  for sample_file in `find $sample_dir -type f -name "*.matrix"`;
  do
    base_name=$(basename $sample_file)
    prefix=${base_name%%.matrix}.${type}
    bsub -P call_m5c -q rhel8_standard -n 11 -M 80000 -W "240:00" \
      -eo $log_folder/$prefix.call_m5c.err -oo $log_folder/$prefix.call_m5c.log \
      -J $prefix.call_m5c \
      "$SCRIPTS/sbatch_call_m5c.sh $sample_file $output_dir $prefix $type $THREAD"
  done
}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

: <<'END'
## for caRNA-batch1
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch1/original_output"
LOG_DIR="$BASE/log/caRNA/call_m5c/batch1"

## exon-based
SAMPLE_DIR="$BASE/sample/caRNA/RNA-m5C/batch1/exon/"
#run_sbatch_call_m5c $SAMPLE_DIR $OUTPUT_DIR $LOG_DIR "exon"

## intron-based
SAMPLE_DIR="$BASE/sample/caRNA/RNA-m5C/batch1/intron/"
#run_sbatch_call_m5c $SAMPLE_DIR $OUTPUT_DIR $LOG_DIR "intron"


## for caRNA-batch2
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch2/original_output"
LOG_DIR="$BASE/log/caRNA/call_m5c/batch2"

create_folder $OUTPUT_DIR
create_folder $LOG_DIR
## exon-based
SAMPLE_DIR="$BASE/sample/caRNA/RNA-m5C/batch2/exon/"
run_sbatch_call_m5c $SAMPLE_DIR $OUTPUT_DIR $LOG_DIR "exon"

## intron-based
SAMPLE_DIR="$BASE/sample/caRNA/RNA-m5C/batch2/intron/"
run_sbatch_call_m5c $SAMPLE_DIR $OUTPUT_DIR $LOG_DIR "intron"
END