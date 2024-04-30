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

def store_feature_stats(feature_dict, feature, repeat, count, peakid_dict, peakid):
    if feature not in feature_dict:
        feature_dict[feature] = defaultdict(float)
    feature_dict[feature][repeat] += count
    return feature_dict

gene_type_list = ["protein_coding", "lncRNA", "pseudogene", "miRNA", "snoRNA", "rRNA", "tRNA", "snRNA", "scRNA", "other_sRNA", 'intergenic']
feature_type_list = ["start_codon", "stop_codon", "5' UTR", "CDS", "3' UTR", "Exon", "Intron", "intergenic"]
#repeat_order_list = ["L1", "L2", "ERVL-MaLR", "ERVL", "ERV1", "Alu", "MIR", "hAT-Charlie", "TcMar-Tigger", "CR1", "Other repeat"]
repeat_order_list = ["Alu", "L1", "MIR", "L2", "ERVL-MaLR", "hAT-Charlie", "ERVL", "TcMar-Tigger", "ERV1", "CR1", "Other repeat"]
order_dict = defaultdict(dict)
count = 0
for i in range(len(gene_type_list)):
    for j in range(len(feature_type_list)):
        gf = '|'.join([gene_type_list[i], feature_type_list[j]])
        order_dict['gf'][gf] = count
        count += 1

for i in range(len(repeat_order_list)):
    order_dict['repeat'][repeat_order_list[i]] = i

name = os.path.splitext(os.path.basename(args.input))[0]

peakid_dict = defaultdict(list)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        peakid = row[3]
        ## record phyloP and length
        if peakid not in peakid_dict:
            peakid_dict[peakid] = defaultdict(list)
        ## record conservation score
        peakid_dict[peakid]['phyloP'] = float(row[4])
        start = int(row[1])
        end = int(row[2])
        peakid_dict[peakid]['length'] = end - start
        ## filter lowly expressed target gene
        geneid = row[7]
        if re.match(r'ENSG', geneid):
            gene_flag =True
        else:
            gene_flag = False
        ##essential information
        gene_type = row[25].split('=')[0]
        tx_type = row[20]
        feature = row[21]
        gf = '|'.join([gene_type, feature])
        peakid_dict[peakid]['gf'].append(gf)
        ##
        if gene_type == "intergenic" and feature == "Exon":
            if tx_type not in order_dict['repeat']:
                repeat = "Other repeat"
            else:
                repeat = tx_type
            peakid_dict[peakid]['repeat'].append(repeat)

feature_dict = defaultdict(dict)
peakid_list = sorted(peakid_dict.keys())
for peakid in peakid_list:
    if 'gf' not in peakid_dict[peakid]:
        feature = 'intergenic'
        count = 1
        feature_dict = store_feature_stats(feature_dict, feature, feature, count, peakid_dict, peakid)
    else:
        gf_list = sorted(peakid_dict[peakid]['gf'], key=lambda x:order_dict['gf'][x])
        if args.mode == "order":
            if gf_list[0] == 'intergenic|intergenic':
                feature = 'intergenic'
                repeat = 'intergenic'
            elif gf_list[0] == 'intergenic|Exon':
                repeat_list = sorted(peakid_dict[peakid]['repeat'], key=lambda x:order_dict['repeat'][x])
                feature = 'intergenic'
                repeat = repeat_list[0]
            elif 'Intron' in gf_list[0]:
                feature = 'Intron'
                if 'repeat' in peakid_dict[peakid]:
                    repeat_list = sorted(peakid_dict[peakid]['repeat'], key=lambda x:order_dict['repeat'][x])
                    repeat = repeat_list[0]
                else:
                    repeat = 'Intron'
            else:
                gene_type, feature = gf_list[0].split('|')
                if gene_type not in ["protein_coding", "IG_gene", "TR_gene"]:
                    if feature == "Exon":
                        feature = "Non-coding exon"
                    else:
                        if feature == "Exon":
                            feature = "Coding exon"
                if feature == "5' UTR":
                    feature = "5-UTR"
                if feature == "3' UTR":
                    feature = "3-UTR"
                repeat = feature
            count = 1
            feature_dict = store_feature_stats(feature_dict, feature, repeat, count, peakid_dict, peakid)
        elif args.mode == "order-average":
            if gf_list[0] == 'intergenic|intergenic':
                feature = 'intergenic'
                repeat = 'intergenic'
                count = 1
                feature_dict = store_feature_stats(feature_dict, feature, repeat, count, peakid_dict, peakid)
            elif gf_list[0] == 'intergenic|Exon':
                feature = 'intergenic'
                repeat_list = sorted(peakid_dict[peakid]['repeat'], key=lambda x:order_dict['repeat'][x])
                for repeat in repeat_list:
                    count = 1 / len(repeat_list)
                    feature_dict = store_feature_stats(feature_dict, feature, repeat, count, peakid_dict, peakid)
            elif 'Intron' in gf_list[0]:
                feature = 'Intron'
                if 'repeat' in peakid_dict[peakid]:
                    repeat_list = ['Intron'] + sorted(peakid_dict[peakid]['repeat'], key=lambda x:order_dict['repeat'][x])
                    for repeat in repeat_list:
                        count = 1 / len(repeat_list)
                        feature_dict = store_feature_stats(feature_dict, feature, repeat, count, peakid_dict, peakid)
                else:
                    count = 1
                    feature_dict = store_feature_stats(feature_dict, feature, feature, count, peakid_dict, peakid)
            else:
                sub_dict = defaultdict(dict)
                for gf in gf_list:
                    gene_type, feature = gf.split('|')
                    if gene_type not in ["protein_coding", "IG_gene", "TR_gene"]:
                        if feature == "Exon":
                            feature = "Non-coding exon"
                    if feature == "5' UTR":
                        feature = "5-UTR"
                    if feature == "3' UTR":
                        feature = "3-UTR"
                    gene_type_order = gene_type_list.index(gene_type)
                    if feature != "intergenic":
                        sub_dict[gene_type_order][feature] = 1
                order = sorted(sub_dict.keys())[0]
                feature_list = sorted(sub_dict[order].keys())
                for feature in feature_list:
                    count = 1 / len(feature_list)
                    feature_dict = store_feature_stats(feature_dict, feature, feature, count, peakid_dict, peakid)
        elif args.mode == "average":
            for gf in gf_list:
                if gene_type not in ["protein_coding", "IG_gene", "TR_gene"]:
                    if feature == "Exon":
                        feature = "Non-coding exon"
                if feature == "5' UTR":
                    feature = "5-UTR"
                if feature == "3' UTR":
                    feature = "3-UTR"
                count = 1 / len(gf_list)
                feature_dict = store_feature_stats(feature_dict, feature, feature, count, peakid_dict, peakid)

feature_type_dict = {"5-UTR":'Exon',
                      "start_codon":'Exon',
                      "CDS":'Exon',
                      "stop_codon":'Exon',
                      "3-UTR":'Exon',
                      "Coding exon":'Exon',
                      "Non-coding exon":'Exon',
                      "Intron":"Intron",
                      "intergenic":"intergenic"}
final_repeat_list = ["L1", "L2", "ERVL-MaLR", "ERVL", "ERV1", "Alu", "MIR", "hAT-Charlie", "TcMar-Tigger", "CR1", "Other repeat"]

with open(args.output, 'w') as out:
    row = ['name', 'type', 'feature', 'repeat', 'count', 'pct\n']
    out.write('\t'.join(row))
    for feature in feature_type_dict:
        main_type = feature_type_dict[feature]
        if feature in feature_dict:
            ## get count and pct for each repeat
            repeat_list = sorted(feature_dict[feature].keys())
            for repeat in repeat_list:
                count = feature_dict[feature][repeat]
                pct = count / len(peakid_list) * 100
                row = [name, main_type, feature, repeat, count, pct]
                out.write('\t'.join(map(str, row)) + '\n')
        else:
            repeat = feature
            count = 0
            pct = 0
            row = [name, main_type, feature, repeat, count, pct]
            out.write('\t'.join(map(str, row)) + '\n')
