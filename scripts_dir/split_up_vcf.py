import argparse
import os

parser = argparse.ArgumentParser(description = 'Split up vcf file into equally sized chunks')
parser.add_argument('vcf_file', type = str, nargs = 1, help = 'vcf file to split up')
parser.add_argument('out_dir', type = str, nargs = 1, help = 'directory to write each new vcf file')
parser.add_argument('prefix', type = str, nargs = 1, help = 'file prefix for each new vcf file')
parser.add_argument('nchunks', type = int, nargs = 1, help = 'number of chunks')
args = parser.parse_args()

if not os.path.isdir(args.out_dir[0]):
    os.mkdir(args.out_dir[0])

with open(args.vcf_file[0], 'r') as data:
    all_file = data.readlines()
    lines_of_header = 0
    for line in all_file:
        lines_of_header += 1
        if line[0:len('#CHROM')] == '#CHROM':
            break
    num_varients = (len(all_file) - lines_of_header)
    lines_per_file = int((num_varients - 1) / args.nchunks[0] + 1)

for i in range(args.nchunks[0]):
    with open(args.out_dir[0] + '/' + args.prefix[0] + '_' + str(i) + '.vcf', 'w') as data:
        data.writelines(all_file[0:lines_of_header])
        for j in range((i * lines_per_file), (min((i + 1) * lines_per_file, num_varients))):
            data.write(all_file[j + lines_of_header])
