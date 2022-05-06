import argparse
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import sys
import os
import json

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
        'chrX'
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
        'chrX' : 156040895
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
args = parser.parse_args()

needed_data = {}
for filenm in args.filenms:
    needed_data[filenm] = []
    with open(filenm, 'r') as data:
        for i, line in enumerate(data):
            if i != 0:
                split_line = line.strip().split('\t')
                if split_line[5] != 'NA':
                    needed_data[filenm].append([split_line[0], int(split_line[1]) + chrom_offset[split_line[0]], int(split_line[2]) + chrom_offset[split_line[0]], float(split_line[5]), split_line[4]])

gene_data = {}
if args.gene_data != '':
    with open(args.gene_data, 'r') as data:
        gene_data = json.load(data)

# make x ticks
xtick_labels = [[chrom, chrom_offset[chrom]] for chrom in chrom_offset]
for gene in gene_data:
    xtick_labels.append([gene, chrom_offset[gene_data[gene][0]] + gene_data[gene][1]])
xtick_labels = sorted(xtick_labels, key = lambda label: label[1])

make_big = 1
fig = plt.figure(figsize = (make_big * 3 * 12, make_big * 2 *12))
axes = fig.subplots(1, 1, sharex = 'all', sharey = 'all')
# axes needs to be a list later in the code so make it one if it isn't already

# count up number of gains and subtract number of losses
gain_or_loss_call = {'NEUT': 0, 'GAIN': 1, 'HETD': -1, 'AMP': 1, 'HLAMP': 1}
bar_values = {}
for j, filenm in enumerate(args.filenms):
    for i, data in enumerate(needed_data[filenm]):
        if data[4] not in gain_or_loss_call:
            sys.exit(f'Erorr: missing {data[4]} in {filenm} from gain_or_loss_call')
        else:
            if j == 0:
                bar_values[f'{data[1]}_{data[2]}'] = gain_or_loss_call[data[4]]
            else:
                bar_values[f'{data[1]}_{data[2]}'] += gain_or_loss_call[data[4]]
bins = [key for key in bar_values]
bins = sorted(bins, key = lambda x: int(x.split('_')[0]))
bar_values = [bar_values[number] for number in bins]

# set colours
colours = []
for value in bar_values:
    if value < 0:
        colours.append('blue')
    elif value > 0:
        colours.append('red')
    elif value ==0:
        colours.append('white')

bars = []
for xs, y, colour in zip(bins, bar_values, colours):
    x1, x2 = xs.split('_')
    x1 = int(x1)
    x2 = int(x2)
    bars.append(Rectangle((x1, 0), x2 - x1, y))

pc = PatchCollection(bars, facecolors = colours)
axes.add_collection(pc)

#min_val = -2 # min([data[2] for data in needed_data[filenm] for filenm in args.filenms[:5]])
#max_val = 2 #max([data[2] for data in needed_data[filenm] for filenm in args.filenms[:5]])
#for i, filenm in enumerate(args.filenms[0:1]):
#    x = [data[1] for data in needed_data[filenm]]
#    y = bar_values
#    axes.scatter(x, y, c = colours, s = 2)
#    axes.set_ylabel(os.path.split(filenm)[1].split('.')[0])
axes.vlines([chrom_offset[chrom] for chrom in chroms], ymin = min(bar_values) - 1, ymax = max(bar_values) + 1)
axes.set_xlim(0, chrom_offset['chrX'] + chr_lengths['chrX'])
axes.set_ylim(bottom = min(bar_values) - 1, top = max(bar_values) + 1)
axes.set_xticks([tick[1] for tick in xtick_labels])
axes.set_xticklabels(tick[0] for tick in xtick_labels)



plt.savefig(args.output)
