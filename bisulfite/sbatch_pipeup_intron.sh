#!/bin/sh

module unload python
module load python/2.7.13
module load hisat/2.1.0
module load samtools/1.17

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/2_m5C_step-by-step_hisat2"
root="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/"

base=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
metadata_dir=$base/metadata

output_dir=$1
prefix=$2

if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

cd $output_dir

echo "pipeup mapped bam"
script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/5_m5C_step-by_step-call_site"

echo "pipeup intron for $prefix.m5C.pileups.tmp"

python $script_dir/m5C_caller_temp_filter.in_memory.mm.intron_v2.py \
  --db $metadata_dir/hg38_gencode_37.noredundance.base \
  --db-intron $metadata_dir/hg38_gencode_37.intron.base \
  -i $prefix.m5C.pileups.tmp \
  -o $prefix.m5C.pileups.formatted.intron.txt \
  --CR $prefix.CR.txt

##
echo "pipeup intron for $prefix.multimappers.m5C.pileups.tmp"

python $script_dir/m5C_caller_temp_filter.in_memory.mm.intron_v2.py \
  --db $metadata_dir/hg38_gencode_37.noredundance.base \
  --db-intron $metadata_dir/hg38_gencode_37.intron.base \
  -i $prefix.multimappers.m5C.pileups.tmp \
  -o $prefix.multimappers.m5C.pileups.formatted.intron.txt \
  --CR $prefix.multimappers.CR.txt
