#!/usr/bin/env python
import os
import sys
import pandas as pd
import math
from glob import glob

folder = sys.argv[1]
output = sys.argv[2]
matrix_list = sorted(glob(os.path.join(folder, '*.matrix')))
con_list = ['WT_rep1', 'WT_rep2', 'WT_rep3', 'KO_rep1', 'KO_rep2', 'KO_rep3', 'KOOE_rep1', 'KOOE_rep2', 'KOOE_rep3']
row_list = []
for matrix in matrix_list:
    batch_name = os.path.basename(matrix).split('.')[0]
    df = pd.read_csv(matrix, sep='\t', header=None, index_col=0)
    for index, dfrow in df.iterrows():
        pval1 = -math.log10(dfrow[10])
        pval2 = -math.log10(dfrow[11])
        for j in range(len(con_list)):
            value = dfrow[j + 1]
            con_name = con_list[j]
            row = [index, batch_name, con_name, value, pval1, pval2]
            row_list.append(row)

output_df = pd.DataFrame(row_list, columns = ['gene', 'batch', 'condition', 'zscore', 'log10pval_KOWT', 'log10pval_OEKO'])
output_df = output_df.set_index('gene')

output_df.to_csv(output, index=True, header=True, sep='\t')
