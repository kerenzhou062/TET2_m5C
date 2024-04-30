#!/usr/bin/env python
import sys
import random

matrix = sys.argv[1]
col = int(sys.argv[2])
maxnum = float(sys.argv[3])
seed = int(sys.argv[4])
random.seed(seed)

with open(matrix, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        row[col] = str(random.uniform(0.01,maxnum))
        print('\t'.join(row))
