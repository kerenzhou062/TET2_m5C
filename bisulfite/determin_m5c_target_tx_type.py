#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='input annotation file')
parser.add_argument('-m', '--mode', action='store', type=str,
                    choices=['order', 'average', 'order-average'],
                    default='order-average',
                    help='how to deal with mulitple annotation of a peak')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def store_feature_stats(feature_dict, feature, count, peakid_dict, peakid):
    if feature not in feature_dict:
        feature_dict[feature] = defaultdict(float)
    feature_dict[feature]['count'] += count
    return feature_dict


gene_type_list = ["snoRNA", "snRNA", "tRNA", "protein_coding", "lncRNA", 'pseudogene', 'intergenic' ]
tx_type_list = ["Alu", "CR1", "ERV1", "ERVL", "ERVK", "ERVL-MaLR", "MIR", "hAT-Charlie", "hAT-Tip100", "hAT-Blackjack", "L1", "L2", "RTE-X", "TcMar-Tigger", "Other repeat", "snRNA", "tRNA"]

name = os.path.splitext(os.path.basename(args.input))[0]

peakid_dict = defaultdict(list)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        peakid = row[3]
        ## filter lowly expressed target gene
        geneid = row[7]
        tx_type = row[20]
        ##essential information
        gene_type = row[25].split('=')[0]
        if gene_type not in gene_type_list:
            peakid_dict[peakid].append('Other')
            continue
        ## skip intergenic
        if tx_type == 'intergenic':
            peakid_dict[peakid].append('Other')
            continue
        ## determin if repeat
        if gene_type == 'intergenic' or tx_type in tx_type_list:
            peakid_dict[peakid].append('repeat')
        else:
            peakid_dict[peakid].append(gene_type)

feature_dict = defaultdict(float)
peakid_list = sorted(peakid_dict.keys())
for peakid in peakid_list:
    each_tx_type_list = peakid_dict[peakid]
    count = 1 / len(each_tx_type_list)
    for tx_type in each_tx_type_list:
        feature_dict[tx_type] += count

all_type_list = sorted(feature_dict.keys())
with open(args.output, 'w') as out:
    row = ['name', 'type', 'count', 'pct\n']
    out.write('\t'.join(row))
    for feature in all_type_list:
        if feature in feature_dict:
            ## get count and pct
            count = feature_dict[feature]
            pct = count / len(peakid_list)
        else:
            count = 0
            pct = 0
        row = [name, feature, count, pct]
        out.write('\t'.join(map(str, row)) + '\n')
