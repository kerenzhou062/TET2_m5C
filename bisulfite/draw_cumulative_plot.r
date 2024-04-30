#!/usr/bin/env Rscript
# exomepeak Script 2
# R script
# Define parameters and load library
library("ggplot2")
library("reshape2")

args = commandArgs(trailingOnly=TRUE)
path <- args[1]
input <- args[2]
pdf <- args[3]
exact <- args[4]
prefix <- args[5]

setwd(path)

CONTROL <- read.table(args[5],sep="\t",header = FALSE)
TREATED <- read.table(args[6],sep="\t",header = FALSE)

pValFunction <- function(pval) {
  if (pval == 0) {
    return ("< 2.22e-308")
  }else if (pval == 1) {
    return ("= 1")
  }else{
    exactVal = sprintf("%.2e", pval)
    exactVal = paste("= ", exactVal, sep = "")
    return (exactVal)
  }
}

myData <- read.table(input,sep="\t", header = TRUE)

if (exact == "TRUE") {
  wilcoxTest <- wilcox.test(TREATED$V1, CONTROL$V1, alternative="less", paired = FALSE, exact = TRUE, correct = TRUE)
}else{
  wilcoxTest <- wilcox.test(TREATED$V1, CONTROL$V1, alternative="less", paired = FALSE, exact = FALSE, correct = FALSE)
}

wilcoxVal = pValFunction(wilcoxTest$p.value)

reshapeData <- melt(myData,id="FC")

levels(reshapeData$variable)[levels(reshapeData$variable)==prefix] <- paste(prefix, " (p ", wilcoxVal,")", sep = "")

teCum <- ggplot(data=reshapeData, aes(x=FC, y=value, colour=variable)) + geom_line(linewidth=1.5)+ labs(x="m5C level", y = "Cumulative frequency") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"), 
    legend.text = element_text(size = 8), 
    legend.title=element_blank(),
    legend.key.width=unit(1,"line"),
    legend.key.height=unit(1,"line"), 
    legend.background=element_rect(fill = "transparent", colour = "transparent"),
    legend.justification=c(0,1), 
    legend.position=c(0.01,0.99), 
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank())

pdf(pdf, height=5)
plot(teCum)
dev.off()
