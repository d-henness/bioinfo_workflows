import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Break up a vcf file by chromosome')
    parser.add_argument('vcf_file', type = str, nargs = 1, help = 'the file name of the vcf file to split')
    parser.add_argument('hla_file', type = str, nargs = 1, help = 'the file name of the hla file')
    parser.add_argument('exp_file', type = str, nargs = 1, help = 'the file name of the expression file')
    args = parser.parse_args()

    path_to_mupexi = '/home/arunimas/MuPeXI/MuPeXI.py'
    path_to_mupexi_config_file = '/home/arunimas/MuPeXI/config.ini'
    with open(args.hla_file[0], 'r') as hla:
        hla_type = hla.readline().split()[1:7]

    with open(args.vcf_file[0], 'r') as data:
        all_lines = data.readlines()
        header = []
        for line in all_lines:
            header.append(line)
            if line[0:4] == '#CHR':
                break

        for i in range(1, 23):
            chromo = 'chr' + str(i)
            filenm = args.vcf_file[0] + '.' + chromo + '.vcf'
            with open(filenm, 'w') as new_vcf:
                for line in header:
                    new_vcf.write(line)
                for line in all_lines:
                    if line.split()[0] == chromo:
                        new_vcf.write(line)
            dirname = args.vcf_file[0] + '.' + chromo
            if not os.path.exists(dirname):
                os.mkdir(dirname)
            base_job_string = """#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=16GB
#SBATCH --account=def-rgni0001

module load miniconda3
source activate mupexi

echo {:}

{:} -v {:} -a {:} -c {:} -t -d {:} -e {:}""".format(filenm, path_to_mupexi, filenm, ','.join(hla_type), path_to_mupexi_config_file, dirname, args.exp_file[0])
            with open(dirname + '.sh', 'w') as job:
                job.write(base_job_string)


        chromo = 'chrX'
        filenm = args.vcf_file[0] + '.' + chromo + '.vcf'
        with open(filenm, 'w') as new_vcf:
            for line in header:
                new_vcf.write(line)
            for line in all_lines:
                if line[0:len(chromo)] == chromo:
                    new_vcf.write(line)
        dirname = args.vcf_file[0] + '.' + chromo
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        base_job_string = """#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=16GB
#SBATCH --account=def-rgni0001

module load miniconda3
source activate mupexi

echo {:}

{:} -v {:} -a {:} -c {:} -t -d {:} -e {:}""".format(filenm, path_to_mupexi, filenm, ','.join(hla_type), path_to_mupexi_config_file, dirname, args.exp_file[0])
        with open(dirname + '.sh', 'w') as job:
            job.write(base_job_string)

        chromo = 'chrY'
        filenm = args.vcf_file[0] + '.' + chromo + '.vcf'
        with open(filenm, 'w') as new_vcf:
            for line in header:
                new_vcf.write(line)
            for line in all_lines:
                if line[0:len(chromo)] == chromo:
                    new_vcf.write(line)
        dirname = args.vcf_file[0] + '.' + chromo
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        base_job_string = """#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=16GB
#SBATCH --account=def-rgni0001

module load miniconda3
source activate mupexi

echo {:}

{:} -v {:} -a {:} -c {:} -t -d {:} -e {:}""".format(filenm, path_to_mupexi, filenm, ','.join(hla_type), path_to_mupexi_config_file, dirname, args.exp_file[0])
        with open(dirname + '.sh', 'w') as job:
            job.write(base_job_string)


if __name__ == '__main__':
    main()
