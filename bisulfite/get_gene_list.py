import sys,os

bed = sys.argv[1]
output = sys.argv[2]
with open(bed, 'r') as f, open(output, 'w') as out:
    for line in f:
        row = line.rstrip().split('\t')
        exon_length = sum([int(i) for i in row[10].strip(',').split(',') if i])
        chrom, start, end =row[0:3]
        strand = row[5]
        start = int(start) - 1
        anno_row = row[3].split(':')
        ##
        gene_id = anno_row[0]
        gene_name = anno_row[1]
        gene_type = anno_row[2]
        tx_id = anno_row[3]
        tx_name = anno_row[4]
        tx_type = anno_row[5]
        output_row = [tx_name, gene_name, tx_id, gene_id, gene_type, chrom, start, end, strand, exon_length]
        out.write('\t'.join(map(str, output_row)) + '\n')
