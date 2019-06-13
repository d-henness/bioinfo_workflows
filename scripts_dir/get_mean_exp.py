import os
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('directory', type = str, help = 'directory containing the bs_*tsv files')
args = parser.parse_args()
if args.directory[-1] != "/":
    args.directory = args.directory + "/"

all_files = os.listdir(args.directory)
bs_files = []
for filenm in all_files:
    if filenm[0:len('bs_')] == 'bs_':
        bs_files.append(filenm)

# read first file to get the number of lines (minus header) per file
with open(args.directory + bs_files[0]) as data:
    lines = data.readlines()
    num_lines = len(lines) - 1
    all_target_ids = []
    for line in lines[1:]:
        all_target_ids.append(line.split()[0])

averages = np.zeros([num_lines, len(bs_files)], dtype=np.float)
for j, filenm in enumerate(bs_files):
    with open(args.directory + filenm) as data:
        header = data.readline() 
        for i, line in enumerate(data):
           averages[i, j] = float(line.split()[-1])

average = averages.mean(1)
print("expression.target_id\tmean_exp\tvar_exp")
for i, var in enumerate(averages.var(1)):
    print(f"{all_target_ids[i]}\t{average[i]}\t{var}")
