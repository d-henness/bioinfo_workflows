import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('map')
parser.add_argument('files', nargs = '+')
args = parser.parse_args()


mapping = {}
with open(args.map, 'r') as data:
    for i, line in enumerate(data):
        if i != 0:
            split_line = line.strip().split()
            mapping[split_line[1].strip('"')] = split_line[2].strip('"')

for filenm in args.files:
    basename = os.path.basename(filenm)
    basename = basename.split('.')
    basename[0] += '-gene_name'
    new_name = '.'.join(basename)

    with open(filenm, 'r') as data:
        with open(new_name, 'w') as new_data:
            for i, line in enumerate(data):
                if i == 0:
                    new_data.writelines(line)
                else:
                    split_line = line.split('\t')
                    if split_line[0] in mapping:
                        split_line[0] = mapping[split_line[0]]
                        new_line = '\t'.join(split_line)
                        new_data.writelines(new_line)
                    else:
                        new_data.writelines(line)
