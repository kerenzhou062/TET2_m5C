#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("ggpubr")

path <- '/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/analysis/qpcr_dotplot'
setwd(path)

spot.theme <- list(
theme_classic(),
theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 15,angle = 90, hjust = 0)),
theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 15)),
theme(legend.position="bottom"), 
theme(axis.line=element_blank()),
theme(text = element_text(size = 15)),
scale_x_discrete(position = "top"))

df <- read.table('THP1_KOOE.combine.txt',sep="\t", header = TRUE)

## re-order batch
df$batch <- factor(df$batch, levels=c('PBS', 'IFN5a', 'IFN5b', 'IFN5r', 'IFN50a', 'IFN50b', 'IFN50r'))
df$condition <- factor(df$condition, levels=c('WT_rep1', 'WT_rep2', 'WT_rep3', 'KO_rep1', 'KO_rep2', 'KO_rep3', 'KOOE_rep1', 'KOOE_rep2', 'KOOE_rep3'))
df$gene <- factor(df$gene, levels=rev(c('OASL', 'ISG15', 'IRF7', 'IFIT2', 'CXCL10')))

plot <- ggplot(df, aes(condition, gene)) + spot.theme +
    geom_point(aes(fill=zscore, size = log10pval_OEKO, colour=log10pval_KOWT), pch=21, stroke=1.5) + 
    scale_fill_gradient2(trans = 'reverse') +
    scale_color_distiller(palette = 'Greens', direction=1) +
    scale_size(range = c(5, 10)) +
    facet_wrap(~ batch)

pdf('THP1_KOOE.combine.dotplot.pdf', width=12, height=10)
plot(plot)
dev.off()

