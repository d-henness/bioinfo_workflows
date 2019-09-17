#!/usr/bin/env python3

import sys
import os
import io
import argparse
import gzip
from zipfile import ZipFile
from operator import itemgetter
import json
import math

def tree_score(tree):
    name, info = tree
    num_samples = len(info['populations']['0']['cellular_prevalence'])
    total_ssms = sum((pop['num_ssms'] for pop in info['populations'].values()))
    normllh_nats = -info['llh'] / total_ssms
    normllh_nats /= num_samples
    normllh_nats /= math.log(2)
    return normllh_nats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--muts', required=True)
    ap.add_argument('-s', '--summ', required=True)
    ap.add_argument('-S', '--ssm_data', required=True)
    ap.add_argument('-t', '--trees', required=True)
    ap.add_argument('vcf', nargs='+')
    ap.add_argument('-o', '--output-dir', required=True)
    args = ap.parse_args()

    muts = json.loads(gzip.open(args.muts, 'rt').read())
    summ = json.loads(gzip.open(args.summ, 'rt').read())

    #put in so ssms can be mapped correctly
    ssms_map = {}
    with open(args.ssm_data, 'r') as data:
        for i, line in enumerate(data):
            if i > 0:
                ssms_map[line.split()[0]] = line.split()[1]


    #tree_list = list(summ['trees'].items())
    #for tup in tree_list:
    #    score = tree_score(tup)
    #    tup[1]['score'] = score
    #sorted_trees = sorted(tree_list, key=lambda x: x[1]['density'])

    tree_densities = [(int(tree), score) for tree, score in summ['tree_densities'].items()]
    sorted_densities = sorted(tree_densities, key=itemgetter(0))  # Sort by tree first so minimum tree idx reported
    sorted_densities = sorted(tree_densities, key=itemgetter(1), reverse=True)
    ssms = muts['ssms']
    print(sorted(ssms.keys()))

    best_tree = str(sorted_densities[0][0])
    best_tree_info = summ['trees'][best_tree]

    with open(f"{args.output_dir}/best_tree.txt", 'w') as data:
        data.write(f"{best_tree}\n")

    with ZipFile(args.trees) as trees_zip:
        best_tree_data = json.loads(trees_zip.open('{}.json'.format(best_tree)).read().decode('utf-8'))
        print(best_tree_data)

    os.makedirs(args.output_dir, exist_ok=True)

    for vcf in args.vcf:
        basename = os.path.basename(vcf)
        lib_dir = f"{args.output_dir}/{basename}"
        os.makedirs(lib_dir, exist_ok=True)

        for pop, assignments in best_tree_data['mut_assignments'].items():
            pop_ssms = assignments['ssms']
            wanted_ssms = set()
            for ssm in pop_ssms:
                #this next line needed to be changed because the needed data is not in the muts file
                chrom, position = ssms_map[ssm].split('_')
                chrom = 'chr' + chrom
                position = int(position)
                wanted_ssms.add((chrom, position))

            population_vcf = "{}/{}.vcf".format(lib_dir, pop)

            with open(vcf, 'r') as vcf_h, open(population_vcf, 'w') as out_h:
                for line in vcf_h:
                    if line.startswith('#'):
                        out_h.write(line)
                        continue
                    sp_line = line.rstrip(os.linesep).split('\t')
                    chrom = sp_line[0]
                    pos = int(sp_line[1])
                    if (chrom, pos) in wanted_ssms:
                        out_h.write(line)


if __name__ == '__main__':
    main()

