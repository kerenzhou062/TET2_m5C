#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("ggpubr")
# or
library("ggsignif")

#args = commandArgs(trailingOnly=TRUE)
#path <- args[1]
#input <- args[2]
#pdf <- args[3]
#prefix <- args[5]

path <- '/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/analysis/violin/batch1'
setwd(path)

data <- read.table('PDX148_EV_TET2.violin.matrix',sep="\t", header = TRUE)

pdf('PDX148_EV_TET2.violin.pdf', width=12, height=8)

plot <- ggplot(data, aes(x = type, y = level, fill = condition)) +
  geom_violin(position = position_dodge(0.8), trim = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.8), outlier.shape = NA)

# Example of adding p-values using ggpubr
plot <- plot + stat_compare_means(method = "wilcox.test", label = "p.signif", aes(group = condition))

plot(plot)
dev.off()

## only L1 + LTR
data <- data[data$type %in% c("L1", "LTR"),]


pdf('PDX148_EV_TET2.L1_LTR.violin.pdf', width=9, height=8)

plot <- ggplot(data, aes(x = type, y = level, fill = condition)) +
  geom_violin(position = position_dodge(0.8), trim = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.8), outlier.shape = NA)

# Example of adding p-values using ggpubr
plot <- plot + stat_compare_means(method = "wilcox.test", label = "p.signif", aes(group = condition))

plot(plot)
dev.off()
