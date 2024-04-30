#!/bin/sh

genome_fa=$1
genome_size=$2
input_bed=$3
output_root_dir=$4
prefix=$5

output_dir=$output_root_dir/$prefix
rm -rf $output_dir
mkdir -p $output_dir

cd $output_dir
base_name=$(basename $input_bed)
base_name=${base_name%.bed}
bedtools slop -i $input_bed -b 5 -g $genome_size -s > ${base_name}.slop.bed
bedtools getfasta -fi $genome_fa -bed ${base_name}.slop.bed -s -rna > ${base_name}.slop.fa

findMotifsGenome.pl \
  $input_bed \
  hg38 \
  $output_dir \
  -rna \
  -noknown \
  -size 20 \
  -len "10" \
  -p 10 \
  -S 15 \
  -minlp "-10" \
  -bits \
  -noweight \
  -preparse \
  -preparsedDir $output_dir/preparse \
  > $output_dir/$prefix.homer.log

rm -rf $output_dir/preparse
