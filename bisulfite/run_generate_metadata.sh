#!/bin/sh
#BSUB -n 1
#BSUB -R "rusage[mem=100GB]"
#BSUB -W 48:00
#BSUB -J run_generate_metadata
#BSUB -P ENCORI
#BSUB -N
#BSUB -oo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/log/run_generate_metadata.log
#BSUB -eo /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/log/run_generate_metadata.err

module unload python
module load python/2.7.13

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/0_m5C_step-by-step_metadata"
root="/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/"
base=$root/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
#hg38
gtf=$root/public/genome/annotation/hg38/v37/gencode.v37.annotation.gtf
bed12=$root/public/genome/annotation/hg38/v37/gencode.v37.annotation.anno.bed12
fasta=$root/public/genome/fasta/hg38/hg38.fa
gsize=$root/public/genome/size/hg38.chrom.sizes
## data
output_dir=$base/metadata

if [[ ! -d "$output_dir" ]]; then
  mkdir -p "$output_dir"
fi

cd  $output_dir
: <<'END'
python $script_dir/gtf2anno.py \
  -i $gtf \
  > hg38_gencode_37.anno

python $base/scripts/get_gene_list.py \
  $bed12 \
  hg38_gencode_37.genelist

python $script_dir/anno_to_base.py \
  -i hg38_gencode_37.anno \
  -o hg38_gencode_37.base

python $script_dir/anno_to_base_remove_redundance_v1.0.py \
  -i hg38_gencode_37.base.sorted \
  -o hg38_gencode_37.noredundance.base \
  -g hg38_gencode_37.genelist
END

script_dir="/research/rgs01/home/clusterHome/kzhou/software/RNA-m5C/5_m5C_step-by_step-call_site"

python $script_dir/anno_to_intronic.py \
  -i hg38_gencode_37.anno \
  -o hg38_gencode_37.intron.base \
  -g hg38_gencode_37.genelist
