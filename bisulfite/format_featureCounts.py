#!/usr/bin/env python3
import os
import sys
import argparse
import re
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input', action='store', type=str, required=True,
                    help='input output from featureCounts')
parser.add_argument('--method', action='store', type=str, choices=['counts', 'FPKM', 'TPM'],
                    default='counts',
                    help='how to estimate gene abundance')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output results')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def normalizeReadCountDf(expMethod, readCountDf):
    lengthDf = readCountDf.iloc[:, 0]
    countDf = readCountDf.iloc[:, range(1, readCountDf.shape[1])]
    totalReadDf = countDf.sum(axis=0)
    if expMethod == 'counts':
        norExpDf = countDf
    elif expMethod == 'FPKM':
        norExpDf = countDf.div(lengthDf, axis=0).div(totalReadDf, axis=1).multiply(1e9)
    elif expMethod == 'TPM':
        countPerbpDf = countDf.div(lengthDf, axis=0)
        norExpDf = countPerbpDf.div(countPerbpDf.sum(axis=0), axis=1).multiply(1e6)
    return norExpDf

##Geneid    Chr Start   End Strand  Length  MMC6_152_input_rep1.genome.sorted.bam   MMC6_152_input_rep2.genome.sorted.bam   MMC6_152_IP_rep1.genome.sorted.bam  MMC6_152_IP_rep2.genome.sorted.bam  MMC6_DMSO_input_rep1.genome.sorted.bam  MMC6_DMSO_input_rep2.genome.sorted.bam  MMC6_DMSO_IP_rep1.genome.sorted.bam MMC6_DMSO_IP_rep2.genome.sorted.bam MOLM13_152_input_rep1.genome.sorted.bam MOLM13_152_input_rep2.genome.sorted.bam MOLM13_152_IP_rep1.genome.sorted.bam    MOLM13_152_IP_rep2.genome.sorted.bam    MOLM13_DMSO_input_rep1.genome.sorted.bam    MOLM13_DMSO_input_rep2.genome.sorted.bam    MOLM13_DMSO_IP_rep1.genome.sorted.bam   MOLM13_DMSO_IP_rep2.genome.sorted.bam
df = pd.read_csv(args.input, sep="\t", header=0, index_col=0, low_memory=False, skiprows=1)
df.drop(columns=['Chr', 'Start', 'End', 'Strand'], axis=1, inplace=True)

## rename the columns
df.rename(columns=lambda x:x.replace(".genome.sorted.bam", ""), inplace=True)

normalizedDf = normalizeReadCountDf(args.method, df)
normalizedDf.to_csv(args.output, header=True, sep="\t")
