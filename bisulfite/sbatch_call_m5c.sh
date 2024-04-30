#!/bin/sh

module unload python
module load python/2.7.13
module load hisat/2.1.0
module load samtools/1.17

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/2_m5C_step-by-step_hisat2"
root="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/"

base=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
hisat2_index=$base/genome/index/hg38_hisat2_index
hisat2_path=/research/rgs01/applications/hpcf/authorized_apps/rhel7_apps/hisat/vendor/2.1.0/
metadata_dir=$base/metadata

sample_file=$1
output_dir=$2
prefix=$3
type=$4
thread=$5

if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

cd $output_dir

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/5_m5C_step-by_step-call_site"

samples=""
while read -r LINE; do
  #dt=( $(perl -p -e 's/\r\n/\n/g' $LINE) )
  LINE=$( echo $LINE | perl -p -e 's/\r\n/\n/g' )
  dt=( $LINE )
  formatted_pipeup=${dt[1]}
  samples="${samples},${formatted_pipeup}"
done < ${sample_file}

if [ "${type}" == "exon" ]; then
  call_script=$script_dir/m5C_caller_multiple.py
else
  call_script=$script_dir/m5C_caller_multiple_v1.4.intron.py
fi

echo "calling peaks on samples: ${samples} "
echo "used calling script: ${call_script} "
python ${call_script} \
  -i $sample_file \
  -o $output_dir/$prefix.tsv \
  -P 2 \
  -c 20 \
  -C 3 \
  -r 0.1 \
  -p 0.05 \
  --processors $thread \
  --method binomial
