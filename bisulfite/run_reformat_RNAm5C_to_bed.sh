#!/bin/sh

function create_folder {
 if [ ! -d "$1" ]; then
   mkdir -p "$1"
  fi
}

m5c_level_cutoff="0.05"

BASE=`realpath ~/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/`
SCRIPTS="$BASE/scripts"
THREAD=10

## for caRNA-batch1
INPUT_DIR="$BASE/m5C_site/caRNA/batch1/original_output"

## for replicates
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch1/bed"
create_folder $OUTPUT_DIR
for m5c_file in `find $INPUT_DIR -type f -name "*.tsv"`;
do
  base_name=$(basename $m5c_file)
  prefix=${base_name%%.tsv}
  $SCRIPTS/reformat_RNAm5C_to_bed.py \
    --input ${m5c_file} \
    --pval 0.05 \
    --type rep \
    --level $m5c_level_cutoff \
    --output $OUTPUT_DIR/${prefix}
done

## merge replicates
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch1/bed"
create_folder $OUTPUT_DIR
for m5c_file in `find $INPUT_DIR -type f -name "*.tsv"`;
do
  base_name=$(basename $m5c_file)
  prefix=${base_name%%.tsv}
  $SCRIPTS/reformat_RNAm5C_to_bed.py \
    --input ${m5c_file} \
    --pval 0.05 \
    --type merge \
    --level $m5c_level_cutoff \
    --output $OUTPUT_DIR/${prefix}.bed
done


## for caRNA-batch2
INPUT_DIR="$BASE/m5C_site/caRNA/batch2/original_output"

## for replicates
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch2/bed"
create_folder $OUTPUT_DIR
for m5c_file in `find $INPUT_DIR -type f -name "*.tsv"`;
do
  base_name=$(basename $m5c_file)
  prefix=${base_name%%.tsv}
  $SCRIPTS/reformat_RNAm5C_to_bed.py \
    --input ${m5c_file} \
    --pval 0.05 \
    --type rep \
    --level $m5c_level_cutoff \
    --output $OUTPUT_DIR/${prefix}
done

## merge replicates
OUTPUT_DIR="$BASE/m5C_site/caRNA/batch2/bed"
create_folder $OUTPUT_DIR
for m5c_file in `find $INPUT_DIR -type f -name "*.tsv"`;
do
  base_name=$(basename $m5c_file)
  prefix=${base_name%%.tsv}
  $SCRIPTS/reformat_RNAm5C_to_bed.py \
    --input ${m5c_file} \
    --pval 0.05 \
    --type merge \
    --level $m5c_level_cutoff \
    --output $OUTPUT_DIR/${prefix}.bed
done
