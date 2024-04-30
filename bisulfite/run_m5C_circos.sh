#!/bin/sh
gsize=/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/public/genome/size/hg38.chrom.sizes
base=/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data
m5c_bed=$base/m5C_site/caRNA/batch1/merge_bed
scripts=$base/scripts
output_dir=$base/analysis/circos/caRNA/PDX148

cd $output_dir
unit=30000000
$scripts/count_5mC_number_density.py \
  --input $m5c_bed/PDX148_EV_bisulfite.merge.bed \
  --size $gsize \
  --unit $unit \
  --output PDX148_EV_bisulfite.merge.circos.txt

$scripts/count_5mC_number_density.py \
  --input $m5c_bed/PDX148_TET2-bisulfite.merge.bed \
  --size $gsize \
  --unit $unit \
  --output PDX148_TET2_bisulfite.merge.circos.txt

min_val=0
max_val=`cat PDX148_EV_bisulfite.merge.circos.txt PDX148_TET2_bisulfite.merge.circos.txt | cut -f 4 | sort -k1,1nr | head -n 1 | awk '{print $1 * 1.01}'`
echo $max_val
## circos path
circos_root=/research/rgs01/applications/hpcf/apps/circos/vendor/0.69

all_ideogram_conf="
<ideogram>\n
\n
<spacing>\n
default = 0.002r\n
break   = 0.2r\n
</spacing>\n
\n
<<include ideogram.position.conf>>\n
<<include ideogram.label.conf>>\n
<<include bands.conf>>\n
\n
radius*       = 0.92r\n
\n
</ideogram>\n

"
all_ideogram_label_conf="
show_label       = yes\n
label_font       = bold\n
label_radius     = dims(image,radius)-95p\n
label_size       = 80\n
label_parallel   = yes\n
label_case       = upper\n
label_format     = eval(sprintf(\"%s\",var(label)))\n
"

config="
<<include $circos_root/etc/colors_fonts_patterns.conf>>\n
\n
<<include ideogram.conf>>\n
<<include ticks.conf>>\n
\n
<image>\n
file* = heatmap.PDX148_EV_TET2_bisulfite.png\n
<<include $circos_root/etc/image.conf>>\n
</image>\n
\n
karyotype   = /research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/metadata/circos/karyotype.human.hg38.txt\n
\n
chromosomes_units  = ${unit}\n
chromosomes        = -hsY\n
#chromosomes_breaks = -hs1:120-140\n
chromosomes_display_default = yes\n
\n
track_width = 0.2\n
track_pad   = 0.02\n
track_start = 0.95\n
\n
<plots>\n
\n
type    = heatmap\n
\n
# default file for all tracks\n
file             = $circos_root/data/karyotype/karyotype.human.txt\n
\n
# a 9 color diverging spectral palette specified using a color list name\n
#color  = oranges-9-seq\n
#color  = grays-9-seq\n
#color  = spectral-9-div-rev
\n
# referenced via conf(plots,color_alt)\n
color_alt = rdylbu-15-div-rev\n
#color_alt = blues-9-seq\n
\n
# or the reverse list\n
#color = greens-6-seq\n
\n
# or you can even combine lists\n
# color = ylorrd-9-seq-rev,ylgnbu-9-seq\n
\n
stroke_thickness = 1\n
stroke_color     = black\n
\n
<plot>\n
  <<include r0r1.conf>>\n
  color            = conf(plots,color_alt)\n
  file             = $output_dir/PDX148_EV_bisulfite.merge.circos.txt\n
  pattern          = hline,solid,vline\n
  color_mapping    = 0  # default\n
  min = ${min_val}\n
  max = ${max_val}\n
  stroke_thickness = 0\n
</plot>\n
\n
<plot>\n
  <<include r0r1.conf>>\n
  color            = conf(plots,color_alt)\n
  file             = $output_dir/PDX148_TET2_bisulfite.merge.circos.txt\n
  pattern          = hline,solid,vline\n
  color_mapping    = 0  # default\n
  min = ${min_val}\n
  max = ${max_val}\n
  stroke_thickness = 0\n
</plot>\n
\n
</plots>\n
\n
<<include $circos_root/etc/housekeeping.conf>>\n

"

module load circos/0.69
my_config_dir=$output_dir/conf/heatmap/
if [ ! -d $my_config_dir ]; then
  mkdir -p $my_config_dir
fi

cd $my_config_dir

echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
echo -e $config > heatmap.PDX148_bisulfite.conf

circos -conf heatmap.PDX148_bisulfite.conf
