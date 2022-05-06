import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('filenm')
parser.add_argument('output')
parser.add_argument('-H', '--height', type = float, default = 0.8)
parser.add_argument('-s', '--sort_order', default = None)
args = parser.parse_args()


all_data = []
samples = []
cell_types = []
with open(args.filenm, 'r') as data:
    for i, line in enumerate(data):
        split_line = line.strip().split('\t')
        if i == 0:
            cell_types = split_line[1:]
        else:
            samples.append(split_line[0])
            all_data.append([float(element) for element in split_line[1:]])

# arrange rows and cols correctly
all_data = np.array(all_data)
all_data = np.transpose(all_data)
sort_order = sorted([i for i, _ in enumerate(cell_types)], key = lambda x: np.average(all_data[x]), reverse = True)
all_data = all_data[sort_order]
cell_types = [cell_types[i] for i in sort_order]

if args.sort_order == None:
    sort_order = sorted([i for i, _ in enumerate(samples)], key = lambda x: all_data[0][x], reverse = True)
else:
    with open(args.sort_order, 'r') as data:
        sample_order = []
        for i, line in enumerate(data):
            sample_order.append(line.split('\t')[0])
    sort_order = [samples.index(sample) for sample in sample_order]
    print(sort_order)

all_data = all_data[:, sort_order]
samples = [samples[i] for i in sort_order]
print(samples)

fig, ax = plt.subplots(figsize = (len(samples) / 4, len(samples) / 4))

norm = colors.Normalize(vmin = 0, vmax = len(cell_types) - 1)
mapper = cm.ScalarMappable(norm = norm, cmap = "PuOr")
left = np.array([0.0 for i in samples])
for i, cell_type in enumerate(cell_types):
    ax.barh(samples, all_data[i], args.height, label = cell_type, left = left, color = mapper.to_rgba(len(cell_types) - i - 1))
    left += all_data[i]
ax.legend()
plt.ylim([-args.height / 2, len(samples) - (args.height / 2)])
plt.savefig(args.output)

plt.clf()
fig, ax = plt.subplots(figsize = (len(samples) / 4, len(samples) / 4))
left = np.array([0.0 for i in samples])
for i, cell_type in enumerate(cell_types):
    if i != 0:
        ax.barh(samples, all_data[i], args.height, label = cell_type, left = left, color = mapper.to_rgba(len(cell_types) - i - 1))
        left += all_data[i]
ax.legend()
plt.ylim([-args.height / 2, len(samples) - (args.height / 2)])
plt.savefig(args.output[:-len(".pdf")] + "_no_other.pdf")
