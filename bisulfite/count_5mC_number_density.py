#!/usr/bin/env python3
import os
import sys
import re
import argparse
import math
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='input bed file')
parser.add_argument('-u', '--unit', action='store', type=int, required=True,
                    help='the unit of density (bp)')
parser.add_argument('-s', '--size', action='store', type=str, required=True,
                    help='genome size file')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

chrom_dict = {}
with open(args.size, 'r') as f:
    for line in f:
        row  = line.strip().split('\t')
        chrom, size = row
        if bool(re.search(r'^chr\d+$', chrom)) or bool(re.search(r'^chr[XYM]$', chrom)):
            chrom_dict[chrom] = [int(int(size) / args.unit) + 1, int(size)]

chrom_unit_count_dict = defaultdict(dict)
with open(args.input, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        unit = int(start / args.unit)
        level = float(row[4])
        if chrom not in chrom_unit_count_dict:
            chrom_unit_count_dict[chrom] = defaultdict(int)
        chrom_unit_count_dict[chrom][unit] += 1

with open(args.output, 'w') as out:
    for chrom in sorted(chrom_dict.keys()):
        unit, size = chrom_dict[chrom]
        for i in range(unit):
            if i not in chrom_unit_count_dict[chrom]:
                count = 0
            else:
                count = chrom_unit_count_dict[chrom][i]
            count = math.log2(count + 1)
            start = i * args.unit
            if i == unit - 1:
                end = size
            else:
                end = (i + 1) * args.unit - 1
            chrom_out = chrom.replace('chr', 'hs')
            row = [chrom_out, start, end, count]
            out.write('\t'.join(map(str, row)) + '\n')
