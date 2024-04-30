#!/bin/sh

input=$1
output=$2
prefix=$3
thread=$4

module load cutadapt/4.4
module load trimmomatic/0.36

cd $output

cutadapt -j $thread --minimum-length 25 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -e 0.25 \
  -q 25 \
  --trim-n \
  -o ${prefix}_R1.cutadapt.fastq.gz -p ${prefix}_R2.cutadapt.fastq.gz \
  $input/${prefix}_R1.fastq.gz $input/${prefix}_R2.fastq.gz

trimmomatic_java=/research/rgs01/applications/hpcf/apps/trimmomatic/vendor/0.36/trimmomatic-0.36.jar

java -jar $trimmomatic_java PE \
  -threads $thread \
  ${prefix}_R1.cutadapt.fastq.gz \
  ${prefix}_R2.cutadapt.fastq.gz \
  ${prefix}.cutadapt.trimmo.rev.fastq.gz \
  ${prefix}.cutadapt.trimmo.revup.fastq.gz \
  ${prefix}.cutadapt.trimmo.fwd.fastq.gz \
  ${prefix}.cutadapt.trimmo.fwdup.fastq.gz \
  HEADCROP:10 SLIDINGWINDOW:4:22 AVGQUAL:25 MINLEN:40
