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

fastq_dir=$1
output_dir=$2
prefix=$3

if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

cd $output_dir

## mapping to the genome
if [ -f "$fastq_dir/${prefix}.cutadapt.trimmo.fwd.fastq" ]; then
  fwd_fastq=$fastq_dir/${prefix}.cutadapt.trimmo.fwd.fastq
else
  cp $fastq_dir/${prefix}.cutadapt.trimmo.fwd.fastq.gz ./
  gunzip ${prefix}.cutadapt.trimmo.fwd.fastq.gz
  fwd_fastq=${prefix}.cutadapt.trimmo.fwd.fastq
fi

if [ -f "$fastq_dir/${prefix}.cutadapt.trimmo.rev.fastq" ]; then
  rev_fastq=$fastq_dir/${prefix}.cutadapt.trimmo.rev.fastq
else
  cp $fastq_dir/${prefix}.cutadapt.trimmo.rev.fastq.gz ./
  gunzip ${prefix}.cutadapt.trimmo.rev.fastq.gz
  rev_fastq=${prefix}.cutadapt.trimmo.rev.fastq
fi

python $script_dir/BS_hisat2.py \
  -F ${fwd_fastq} \
  -R ${rev_fastq} \
  -o $prefix \
  -I $hisat2_index \
  --index-prefix HISAT2 \
  --hisat2-path $hisat2_path \
  --del-convert \
  --del-sam

## delete fastq
rm -f *.fastq

## sort and index mapped bam
echo "sorting and indexing $prefix.bam"
samtools sort $prefix.bam -o $prefix.sorted.bam --threads 15
samtools index $prefix.sorted.bam

echo "sorting and indexing $prefix.multimappers.bam"
samtools sort $prefix.multimappers.bam -o $prefix.multimappers.sorted.bam --threads 15
samtools index $prefix.multimappers.sorted.bam

echo "pipeup mapped bam"
script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/4_m5C_step-by-step_pileup"
fasta=$root/public/genome/fasta/hg38/hg38.fa

echo "pipeup $prefix.sorted.bam"
python $script_dir/pileup_genome_multiprocessing_v1.4.py \
  -P 10 \
  -f $fasta \
  -i $prefix.sorted.bam \
  -o $prefix.m5C.pileups.tmp

python $script_dir/m5C_pileup_formatter.py \
  --db $metadata_dir/hg38_gencode_37.noredundance.base \
  -i $prefix.m5C.pileups.tmp \
  -o $prefix.m5C.pileups.formatted.txt \
  --CR $prefix.CR.txt

##
echo "pipeup $prefix.multimappers.sorted.bam"
python $script_dir/pileup_genome_multiprocessing_v1.4.py \
  -P 10 \
  -f $fasta \
  -i $prefix.multimappers.sorted.bam \
  -o $prefix.multimappers.m5C.pileups.tmp

python $script_dir/m5C_pileup_formatter.py \
  --db $metadata_dir/hg38_gencode_37.noredundance.base \
  -i $prefix.multimappers.m5C.pileups.tmp \
  -o $prefix.multimappers.m5C.pileups.formatted.txt \
  --CR $prefix.multimappers.CR.txt

echo "pipeup mapped bam for intron"
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
