#!/bin/sh
module load meme/5.5.4

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
bedtools slop -i $input_bed -b 10 -g $genome_size -s > ${base_name}.slop.bed
bedtools getfasta -fi $genome_fa -bed ${base_name}.slop.bed -s -rna | awk '{seq=toupper($0);seq=gensub(/T/, "U", "g", seq);print seq}'> ${base_name}.slop.fa

meme -rna -objfun ce -nmotifs 5 -cefrac 0.1 -hsfrac 0.7 -minw 8 -maxw 10 -seed 100 -oc ./ ${base_name}.slop.fa
