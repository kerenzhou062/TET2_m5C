#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import scipy.stats as stats
import pandas as pd
import numpy as np
import math

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

def store_rna_location_stats(rna_location_dict, location_list, count, peakid_dict, peakid):
    location_count = len(location_list)
    each_count = count / location_count
    for location in location_list:
        if location not in rna_location_dict:
            rna_location_dict[location] = defaultdict(float)
        rna_location_dict[location]['count'] += each_count
        rna_location_dict[location]['phyloP'] += peakid_dict[peakid]['phyloP'] * peakid_dict[peakid]['length']
        rna_location_dict[location]['length'] += peakid_dict[peakid]['length']
    return rna_location_dict

name = os.path.splitext(os.path.basename(args.input))[0]

final_rna_location_list = []
peakid_dict = defaultdict(dict)
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
        gene_type = row[25]
        rna_sublocation = row[26]
        ## skip intergenic
        if gene_type == 'intergenic':
            continue
        rna_sublocation_list = rna_sublocation.split(',')
        peakid_dict[peakid]['location'].append(rna_sublocation_list)
        final_rna_location_list.extend(rna_sublocation_list)

##
rna_location_dict = defaultdict(dict)
peakid_list = sorted(peakid_dict.keys())
for peakid in peakid_list:
    if 'location' not in peakid_dict[peakid]:
        location_list = ['Unkown']
        count = 1
        rna_location_dict = store_rna_location_stats(rna_location_dict, location_list, count, peakid_dict, peakid)
    else:
        rna_location_list_count = len(peakid_dict[peakid]['location'])
        count = 1 / rna_location_list_count
        for location_list in peakid_dict[peakid]['location']:
            rna_location_dict = store_rna_location_stats(rna_location_dict, location_list, count, peakid_dict, peakid)

final_rna_location_list = sorted(set(final_rna_location_list))

with open(args.output, 'w') as out:
    row = ['rbp_label', 'type', 'count', 'pct', 'phyloP', 'phyloP_height\n']
    out.write('\t'.join(row))
    for location in final_rna_location_list:
        if location in rna_location_dict:
            ## get count and pct
            count = rna_location_dict[location]['count']
            pct = count / len(peakid_list)
            ## get average conservation score
            phyloP = rna_location_dict[location]['phyloP'] / rna_location_dict[location]['length']
            phyloP_height = 1 / len(final_rna_location_list)
        else:
            count = 0
            pct = 0
            phyloP = 0
            phyloP_height = 1 / len(final_rna_location_list)
        row = [name, location, count, pct, phyloP, phyloP_height]
        out.write('\t'.join(map(str, row)) + '\n')
