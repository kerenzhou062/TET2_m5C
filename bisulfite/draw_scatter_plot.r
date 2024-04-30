#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")
library('ggrepel')

#args = commandArgs(trailingOnly=TRUE)
#path <- args[1]
#input <- args[2]
#pdf <- args[3]
#prefix <- args[5]

path <- '/research_jude/rgs01_jude/groups/tatevgrp/home/kzhou/project/chen_yang_lab/ying_qing/tet1/bisulfite_seq/inhouse_data/analysis/scatter/batch1'
setwd(path)

myData <- read.table('PDX148_EV_TET2.scatter.label.matrix',sep="\t", header = TRUE)

myData$size <- 4
myData[ (myData$x < 0.05) | (myData$y < 0.05), ]$size = 4

myData$alpha <- 0.7
myData[ (myData$x < 0.05) | (myData$y < 0.05), ]$alpha = 0.7

subData <- myData[ ! is.na(myData$gene), ]

pdf('PDX148_EV_TET2.scatter.label.pdf', width=9, height=8)

plot <- ggplot(myData, 
    aes(x = x, y = y, color = status)) +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_abline(slope = 1,
            intercept = 0,
            color="gray",
            linetype = "dashed") +
    geom_point(size=4, alpha=0.7) +
    geom_point(data=subData, size=4, alpha=0.8, color='#1985A1') +
    scale_color_manual(values = c("Hypo" = "#66C2A5", "Hyper" = "#FC8D62", "nochange"="gray")) +
    xlab("m5C level in PDX148 EV cells") + 
    ylab("m5C level in PDX148 TET2-OE cells") +
    geom_label_repel(data=subset(myData, gene != "NA"),
                  aes(label=gene),
                  box.padding = 0.25,
                  label.padding = 0.25,
                  min.segment.length = 0,
                  force_pull=1,
                  segment.color = 'grey50',
                  seed = 100,
                  show.legend=FALSE) +
    theme_bw() +
    theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
    )

plot(plot)
dev.off()
