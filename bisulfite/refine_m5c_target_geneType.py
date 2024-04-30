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
parser.add_argument('-a', '--anno', action='store', type=str, required=True,
                    help='more detail annotation file  include repeats (gene_type.txt)')
parser.add_argument('-p', '--protein', action='store', type=str, required=True,
                    help='classified proteins matrix collected from protein-atlas')
parser.add_argument('-r', '--rnasub', action='store', type=str, required=True,
                    help='All experimental RNA subcellular localization data from RNALocate v2')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

anno_dict = {}
with open(args.anno, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        gene_type = row[0]
        sub_class = row[1]
        main_class = row[2]
        anno_dict[gene_type] = [sub_class, main_class]

protein_dict = defaultdict(list)
with open(args.protein, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        ## no version
        geneid = row[0]
        protein_class = row[1]
        if protein_class in ['TCR', 'Immunoglobulin']:
            protein_class = "Other_protein"
        protein_dict[geneid].append(protein_class)
## remove duplicates
for geneid in sorted(protein_dict.keys()):
    protein_dict[geneid] = sorted(set(protein_dict[geneid]))

rna_locate_dict = defaultdict(list)
with open(args.rnasub, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        ##Homo sapiens
        species = row[5]
        if species != "Homo sapiens":
            continue
        gene_type = row[4]
        if gene_type == "circRNA":
            continue
        ## symbol
        gene_name = row[3]
        location = row[7]
        rna_locate_dict[gene_name].append(location)

## remove duplicates
rna_locate_stats_dict = defaultdict(int)
for gene_name in sorted(rna_locate_dict.keys()):
    rna_locate_list = sorted(set(rna_locate_dict[gene_name]))
    rna_locate_dict[gene_name] = rna_locate_list
    for rna_locate in rna_locate_list:
        rna_locate_stats_dict[rna_locate] += 1

## only pick the top 20
rna_locate_selected_list = sorted(rna_locate_stats_dict.keys(), key=lambda x:-rna_locate_stats_dict[x])[0:20]
## re name the rna_location
for gene_name in sorted(rna_locate_dict.keys()):
    rna_locate_list = rna_locate_dict[gene_name]
    new_locate_list = []
    for rna_locate in rna_locate_list:
        if rna_locate not in rna_locate_selected_list:
            new_locate_list.append("Other")
        else:
            new_locate_list.append(rna_locate)
    rna_locate_dict[gene_name] = sorted(set(new_locate_list))

## re annote
with open(args.output, 'w') as out, open(args.input, 'r') as f:
    row = f.readline().strip().split('\t')
    row += ['GeneTypeDetail', "RnaLocation"]
    out.write('\t'.join(row) + '\n')
    for line in f:
        row = line.strip().split('\t')
        gene_type = row[16]
        class_type = row[17]
        tx_type = row[20]
        ## remove version
        geneid = row[12].split('.')[0]
        ## determin detail of gene typ
        if class_type == 'other' and bool(re.match(r'ENSG', geneid)) is False:
            detail_type, main_type = anno_dict[gene_type]
            if detail_type in ['pRNA', 'paRNA', 'eRNA']:
                detail_type = 'intergenic'
                main_type = 'intergenic'
        else:
            detail_type, main_type = anno_dict[gene_type]
        ## if gene in protein class, then re label as protein_coding
        if geneid in protein_dict:
            detail_type = 'protein_coding'
        ## mark repeat as intergenic
        if main_type == 'repeat':
            detail_type = 'intergenic'
        ## if protein_coding
        if detail_type == 'protein_coding':
            if geneid in protein_dict:
                protein_class = ','.join(protein_dict[geneid])
            else:
                protein_class = 'Other_protein'
            detail_type = 'protein_coding=' + protein_class
        ## rna subcellular location
        gene_name = row[13]
        if gene_name in rna_locate_dict:
            rna_location = ','.join(rna_locate_dict[gene_name])
        else:
            rna_location = 'Unkown'
        row += [detail_type, rna_location]
        out.write('\t'.join(row) + '\n')
