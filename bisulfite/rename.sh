#!/bin/sh

function matrixToRun {
  ## $PARAM1 and $PARAM1 used for filter
  SAMPLE_MATRIX=$1
  i=1
  while IFS=$'\t' read -r -a matrixArr
  do
    test $i -eq 1 && ((i=i+1)) && continue
    original_id="${matrixArr[1]}"
    exp_prefix="${matrixArr[0]}"
    rename "$original_id" "${exp_prefix}" ${original_id}*.fastq.gz
  done < $SAMPLE_MATRIX
}

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`

: <<'END'
## for caRNA-batch1
FASTQ_DIR="$BASE/fastq/caRNA/batch1/original"
cd $FASTQ_DIR
matrixToRun $BASE/sample/caRNA/sample.batch1.matrix
END

## for caRNA-batch1
FASTQ_DIR="$BASE/fastq/caRNA/batch2/original"
cd $FASTQ_DIR
matrixToRun $BASE/sample/caRNA/sample.batch2.matrix
