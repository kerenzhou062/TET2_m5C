#!/usr/bin/env python3
import os
import sys
import re
import argparse
import math
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', type=str, required=True,
                    help='input RBP-table matrix')
parser.add_argument('--type', action='store', type=str,
                    choices=['merge', 'rep'],
                    help='whether output sites that only significant in replicates')
parser.add_argument('--level', action='store', type=float,
                    default=0.1,
                    help='cutoff of m5C level for m5C site')
parser.add_argument('--pval', action='store', type=float,
                    default=0.05,
                    help='cutoff for m5C site')
parser.add_argument('--output', action='store', type=str, required=True,
                    help='output bed')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

'''
                  0                   1         2     3    4         5             6    7        8     9    10        11            12   13       14
0  chr10@100356731@+   ENSG00000099194.6       SCD  22.0  1.0  0.045455  2.070159e-02  1.0  0.00095  20.0  3.0  0.150000  3.296035e-06  1.0  0.00143
1  chr10@103178170@-                 NaN       NaN  24.0  3.0  0.125000  1.085822e-05  1.0  0.00177  15.0  0.0  0.000000  1.000000e+00  1.0  0.00158
2  chr10@103801300@-  ENSG00000107957.16  SH3PXD2A  24.0  3.0  0.125000  3.317070e-06  1.0  0.00119  24.0  0.0  0.000000  1.000000e+00  1.0  0.00123
3  chr10@103838106@-                 NaN       NaN  14.0  0.0  0.000000  1.000000e+00  1.0  0.00177  26.0  4.0  0.153846  9.072386e-08  1.0  0.00158
4  chr10@109882554@-  ENSG00000108039.18   XPNPEP1  28.0  4.0  0.142857  8.211539e-08  1.0  0.00142  25.0  1.0  0.040000  3.279980e-02  1.0  0.00133
'''
df = pd.read_csv(args.input, sep='\t', skiprows=2, header=None)
if args.type == 'merge':
    df = df[(df.iloc[:,6] < args.pval) & (df.iloc[:,12] < args.pval)]
    ## get average of coverage
    df[15] = df.iloc[:,[3, 9]].mean(axis=1)
    ## get average C count
    df[16] = df.iloc[:,[4, 10]].mean(axis=1)
    ## get m5c level
    df[17] = df[16].divide(df[15])

    ## reformat m5C site
    m5C_list = list(map(lambda x:x.split('@'), df.iloc[:,0].to_list()))
    chrom_list = list(map(lambda x:x[0], m5C_list))
    start_list = list(map(lambda x:int(x[1]) - 1, m5C_list))
    end_list = list(map(lambda x:x[1], m5C_list))
    strand_list = list(map(lambda x:x[2], m5C_list))

    start_list = [ 0 if x < 0 else x for x in start_list]

    df[18] = chrom_list
    df[19] = start_list
    df[20] = end_list
    df[21] = strand_list

    final_df = df.loc[:,[18, 19, 20, 1, 17, 21, 2, 15, 16, 6, 12]]
    final_df.columns = ['#chrom', 'start', 'end', 'gene_id', 'm5c_level', 'strand', 'gane_name', 'total_reads_count', 'C_reads_count', 'pval_rep1', 'pval_rep2']
    final_df = final_df.fillna('NA')
    ## filter low m5c_level
    final_df = final_df.loc[(final_df['m5c_level'] >= args.level)]
    final_df.to_csv(args.output, sep='\t', index=False, header=True)
else:
    index_list = [[0,1,2,3,4,5,6,7,8], [0,1,2,9,10,11,12,13,14]]
    for i in range(2):
        sub_df = df.loc[:,index_list[i]]
        sub_df.columns = index_list[0]
        sub_df = sub_df[(sub_df.iloc[:,6] < args.pval)]
        ## reformat m5C site
        m5C_list = list(map(lambda x:x.split('@'), sub_df.iloc[:,0].to_list()))
        chrom_list = list(map(lambda x:x[0], m5C_list))
        start_list = list(map(lambda x:int(x[1]) - 1, m5C_list))
        end_list = list(map(lambda x:x[1], m5C_list))
        strand_list = list(map(lambda x:x[2], m5C_list))
        start_list = [ 0 if x < 0 else x for x in start_list]
        sub_df[9] = chrom_list
        sub_df[10] = start_list
        sub_df[11] = end_list
        sub_df[12] = strand_list
        final_df = sub_df.loc[:,[9, 10, 11, 1, 5, 12, 2, 3, 4, 6]]
        final_df.columns = ['#chrom', 'start', 'end', 'gene_id', 'm5c_level', 'strand', 'gane_name', 'total_reads_count', 'C_reads_count', 'pval']
        ## filter low m5c_level
        final_df = final_df[final_df['m5c_level'].notnull()]
        final_df = final_df.fillna('NA')
        final_df = final_df.loc[(final_df['m5c_level']>= args.level)]
        ## output bed per replicate
        rep = 'rep{}'.format(i + 1)
        output = '{}.{}.bed'.format(args.output, rep)
        final_df.to_csv(output, sep='\t', index=False, header=True)
