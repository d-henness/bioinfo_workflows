import argparse
import matplotlib.pyplot as plt
import re
import sys
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('catagory')
parser.add_argument('-f', '--font_size', type = int, default = 2)
args = parser.parse_args()

starting_signature_index = None

with open(args.input, 'r') as data:
    fields = []
    signatures = []
    sample_data = []
    for i, line in enumerate(data):
        if i == 0:
            fields = line.strip().split('\t')
            for j, field in enumerate(fields):
                print(field)
                if field.startswith('S1-'):
                    signatures = fields[j:]
                    starting_signature_index = j
                    break
        else:
            split_line = line.strip().split('\t')
            for i in range(starting_signature_index, len(split_line)):
                split_line[i] = float(split_line[i])
            sample_data.append(split_line)


    print(len(sample_data))
    sum_of_signatures_per_sample = {sample[0]: sum(signature for signature in sample[starting_signature_index:]) for i, sample in enumerate(sample_data)}
    print(len(sum_of_signatures_per_sample))

    sample_data = sorted(sample_data, key = lambda x: (x[fields.index(args.catagory)], sum_of_signatures_per_sample[x[0]] ), reverse = True)
    sample_names = [sample[0] for sample in sample_data]
    divisions = []
    division_labels = []
    group = None
    for i, sample in enumerate(sample_data):
        if sample[fields.index(args.catagory)] != group:
            group = sample[fields.index(args.catagory)]
            divisions.append(i)
            division_labels.append(group)

fig, ax = plt.subplots(1,2)

for i, signal in enumerate(signatures):
    ax[0].barh(
        [j for j, _ in enumerate(sample_data)],
        [sample_data[j][i + starting_signature_index] for j, _ in enumerate(sample_names)],
        left = [sum(sample_data[j][starting_signature_index: i + starting_signature_index]) for j, _ in enumerate(sample_names)],
        label = signal
    )
ax[0].set_xlabel('Counts')
ax[0].set_yticks(divisions)
ax[0].set_yticklabels(division_labels, fontsize = args.font_size)
ax[0].legend()


for i, signal in enumerate(signatures):
    data = [sample_data[j][i + starting_signature_index] / sum_of_signatures_per_sample[sample] for j, sample in enumerate(sample_names)]
    lefts = [sum(sample_data[j][starting_signature_index: i + starting_signature_index]) / sum_of_signatures_per_sample[sample] for j, sample in enumerate(sample_names)]
    print(len(data))
    ax[1].barh(
        [j for j, _ in enumerate(sample_data)],
        data,
        left = lefts,
        label = signal
    )
ax[1].set_xlabel('Fraction')
ax[1].set_yticks(divisions)
ax[1].set_yticklabels(division_labels, fontsize = args.font_size)

plt.savefig(f'{args.output}.pdf')
