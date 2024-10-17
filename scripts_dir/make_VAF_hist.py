#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
import os

def main():
    filenm = sys.argv[1]
    outfilenm = f'{os.path.splitext(filenm)[0]}_hist.pdf'

    all_vafs = []
    with open(filenm, 'r') as data:
        for i, line in enumerate(data):
            split_line = line.strip().split('\t')
            vafs = split_line[-1].split(',')
            for vaf in vafs:
                all_vafs.append(float(vaf))

    n, bins, patches = plt.hist(all_vafs, bins = 50)
    plt.xlabel('VAF')
    plt.ylabel('Counts')
    plt.ylim(bottom = 1)
    plt.savefig(outfilenm)



if __name__ == "__main__":
    main()
