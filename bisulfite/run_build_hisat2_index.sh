#!/bin/sh
#BSUB -n 21
#BSUB -R "rusage[mem=15GB]"
#BSUB -W 48:00
#BSUB -J run_build_hisat2_index
#BSUB -P ENCORI
#BSUB -N
#BSUB -oo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/log/run_build_hisat2_index.log
#BSUB -eo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/log/run_build_hisat2_index.err

module unload python
module load python/2.7.13
module load hisat/2.1.0

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/0_m5C_step-by-step_metadata"
root="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/"
base=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
#hg38
gtf=$root/public/genome/annotation/hg38/v37/gencode.v37.annotation.gtf
fasta=$root/public/genome/fasta/hg38/hg38.fa
gsize=$root/public/genome/size/hg38.chrom.sizes
## data
OUTPUT_DIR=$base/genome/index
HISAT2_PATH=/research/rgs01/applications/hpcf/authorized_apps/rhel7_apps/hisat/vendor/2.1.0/

python $script_dir/BS_hisat2_index.py \
  --threads 20 \
  --input $fasta \
  --gtf $gtf \
  --hisat2-path $HISAT2_PATH \
  --output $OUTPUT_DIR/hg38_hisat2_index
