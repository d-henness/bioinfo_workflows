import argparse
import matplotlib.pyplot as plt
import sys
import os
import json
import collections

chroms = [
        'chr1',
        'chr2',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr7',
        'chr8',
        'chr9',
        'chr10',
        'chr11',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        'chr17',
        'chr18',
        'chr19',
        'chr20',
        'chr21',
        'chr22',
        'chrX',
        'total'
]

chr_lengths = {
        'chr1':  248956422,
        'chr2':  242193529,
        'chr3':  198295559,
        'chr4':  190214555,
        'chr5':  181538259,
        'chr6':  170805979,
        'chr7':  159345973,
        'chr8':  145138636,
        'chr9':  138394717,
        'chr10': 133797422,
        'chr11': 135086622,
        'chr12': 133275309,
        'chr13': 114364328,
        'chr14': 107043718,
        'chr15': 101991189,
        'chr16':  90338345,
        'chr17':  83257441,
        'chr18':  80373285,
        'chr19':  58617616,
        'chr20':  64444167,
        'chr21':  46709983,
        'chr22':  50818468,
        'chrX' : 156040895,
        'total': 0
}

chrom_offset = {}
for chrom in chroms:
    if chrom == 'chr1':
        chrom_offset[chrom] = 0
    else:
        chrom_offset[chrom] = sum([chr_lengths[chro] for chro in chroms[0:chroms.index(chrom)]])


parser = argparse.ArgumentParser()
parser.add_argument('filenms', nargs = '+')
parser.add_argument('-o', '--output', default = 'test.pdf')
parser.add_argument('-g', '--gene_data', default = '')
parser.add_argument('-b', '--bin_width', type = int, default = 1000000)
args = parser.parse_args()


#make the bin widths

bin_widths = []
i = 0
for chrom in chroms:
    while i < chrom_offset[chrom]:
        bin_widths.append(i)
        i += args.bin_width
    i = chrom_offset[chrom]

# read the data from the files
needed_data = []
missense = []
for filenm in args.filenms:
    with open(filenm, 'r') as data:
        for i, line in enumerate(data):
            if i != 0:
                split_line = line.strip().split('\t')
                if split_line[0] == 'chrY':
                    split_line[0] = 'chrX'
                needed_data.append(chrom_offset[split_line[0]] + int(split_line[1]))
                if 'missense_variant' in line:
                    missense.append(chrom_offset[split_line[0]] + int(split_line[1]))


gene_data = {}
if args.gene_data != '':
    with open(args.gene_data, 'r') as data:
        gene_data = json.load(data)

# make x ticks
xtick_labels = [[chrom, chrom_offset[chrom]] for chrom in chroms[:-1]]
for gene in gene_data:
    xtick_labels.append([gene, chrom_offset[gene_data[gene][0]] + gene_data[gene][1]])
xtick_labels = sorted(xtick_labels, key = lambda label: label[1])

make_big = 1
fig = plt.figure(figsize = (make_big * 3 * 12, make_big * 2 * 3))
axes = fig.subplots(1, 1, sharex = 'all', sharey = 'all')
# axes needs to be a list later in the code so make it one if it isn't already
#if len(args.filenms) == 1:
#    axes = [axes]

#for i, filenm in enumerate(args.filenms):
values, _, _ = axes.hist(needed_data, bin_widths, label = 'All Mutations')
axes.hist(missense, bin_widths, label = 'Missense Mutations')
axes.set_ylabel('Mutations')
axes.vlines([chrom_offset[chrom] for chrom in chroms], 0, max(values) * 1.1)
axes.set_xlim(0, chrom_offset['total'])
axes.set_ylim(bottom = 0, top = max(values) * 1.1)
axes.set_xticks([tick[1] for tick in xtick_labels])
axes.set_xticklabels(tick[0] for tick in xtick_labels)



plt.legend()
plt.savefig(args.output)
