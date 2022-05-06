import sys

gene_ids = sys.argv[1]
filenm = sys.argv[2]

gene_map = {}
with open (gene_ids, 'r') as data:
    for line in data:
        split_line = line.strip().split()
        gene_map[split_line[0]] = split_line[1]

with open(filenm, 'r') as data:
    for line in data:
        split_line = line.strip().split('\t')
        if split_line[0] in gene_map:
            split_line[0] = gene_map[split_line[0]]
            new_line = '\t'.join(split_line)
            print(new_line)
        else:
            print(line)
