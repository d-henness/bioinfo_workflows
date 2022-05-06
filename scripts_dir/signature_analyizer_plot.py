import argparse
import matplotlib.pyplot as plt
import re
import sys
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('output')
parser.add_argument('-f', '--font_size', type = int, default = 16)
args = parser.parse_args()

sample_name_index = None
signature_indices = []
signature_names = []

with open(args.input, 'r') as data:
    valid = re.compile(r"^S[0-9]+-")
    split_line = data.readline().split('\t')
    for i, element in enumerate(split_line):
        if element == "Tumor_Sample_Barcode":
            sample_name_index = i
        elif re.match(valid, element) != None:
            signature_indices.append(i)
            signature_names.append(element)

    if sample_name_index == None:
        sys.exit("Could not find Tumor_Sample_Barcode column")
    elif len(signature_names) == 0:
        sys.exit("Could not find signatures")

    print(signature_indices)

    sum_of_signatures_per_sample = defaultdict(lambda: [0 for i in signature_indices])

    for line in data:
        split_line = line.split('\t')
        for i, index in enumerate(signature_indices):
            try:
                sum_of_signatures_per_sample[split_line[sample_name_index]][i] += float(split_line[index])
            except:
                continue

    print(sum_of_signatures_per_sample)

    sample_names = list(sum_of_signatures_per_sample.keys())
    sample_names = sorted(sample_names, key = lambda x: sum(sum_of_signatures_per_sample[x]), reverse = True)
    sample_names_parsed = [re.sub(r'^RG', r'MF', sample_name) for sample_name in sample_names]
    sample_names_parsed = [re.sub(r'-', r'_', sample_name) for sample_name in sample_names_parsed]
    sample_names_parsed = [re.sub(r'_0$', r'', sample_name) for sample_name in sample_names_parsed]
    print(sample_names)

with open(f'{args.output}.tsv', 'w') as sum_data:
    sig_string = '\t'.join(signature_names)
    sum_data.write(f'sample\t{sig_string}')
    for i, sample in enumerate(sample_names):
        line_string = '\t'.join([str(i) for i in sum_of_signatures_per_sample[sample]])
        sum_data.write(f'{sample_names_parsed[i]}\t{line_string}\n')



    fig, ax = plt.subplots(1,2)

    for i, signal in enumerate(signature_names):
        ax[0].barh(
            sample_names,
            [sum_of_signatures_per_sample[sample][i] for sample in sample_names],
            left = [sum(sum_of_signatures_per_sample[sample][:i]) for sample in sample_names],
            label = signal
        )
    ax[0].set_ylabel('Counts')
    ax[0].set_yticklabels(labels = sample_names_parsed, fontsize = args.font_size)
    ax[0].legend()

    muts_per_sample = [sum(sum_of_signatures_per_sample[sample]) for sample in sample_names]
    for i, signal in enumerate(signature_names):
        ax[1].barh(
            sample_names,
            [(sum_of_signatures_per_sample[sample][i] / muts_per_sample[j]) if muts_per_sample[j] != 0 else 0 for j, sample in enumerate(sample_names)],
            left = [(sum(sum_of_signatures_per_sample[sample][:i]) / muts_per_sample[j]) if muts_per_sample[j] != 0 else 0 for j, sample in enumerate(sample_names)],
            label = signal
        )
    ax[1].set_ylabel('Fraction')
    ax[1].set_yticklabels(labels = sample_names_parsed, fontsize = args.font_size)

    plt.savefig(f'{args.output}.pdf')
