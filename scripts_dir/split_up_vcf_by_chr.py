import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Break up a vcf file by chromosome')
    parser.add_argument('vcf_file', type = str, nargs = 1, help = 'the file name of the vcf file to split')
    parser.add_argument('out_dir', type = str, help = 'the directory to write vcf files to')
    parser.add_argument('prefix', type = str, help = 'prefix of vcf file')
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    with open(args.vcf_file[0], 'r') as data:
        all_lines = data.readlines()
        header = []
        for line in all_lines:
            header.append(line)
            if line[:len('#CHR')] == '#CHR':
                break

        for i in range(1, 23):
            chromo = f'chr{i}'
            filenm = f"{args.out_dir}/{args.prefix}_{chromo}.vcf"
            if not os.path.exists(f"{args.out_dir}/{chromo}"):
                os.mkdir(f"{args.out_dir}/{chromo}")
            with open(filenm, 'w') as new_vcf:
                for line in header:
                    new_vcf.write(line)
                for line in all_lines:
                    if line.split()[0] == chromo:
                        new_vcf.write(line)

        chromo = 'chrX'
        filenm = f"{args.out_dir}/{args.prefix}_{chromo}.vcf"
        if not os.path.exists(f"{args.out_dir}/{chromo}"):
            os.mkdir(f"{args.out_dir}/{chromo}")
        with open(filenm, 'w') as new_vcf:
            for line in header:
                new_vcf.write(line)
            for line in all_lines:
                if line.split()[0] == chromo:
                    new_vcf.write(line)

        chromo = 'chrY'
        filenm = f"{args.out_dir}/{args.prefix}_{chromo}.vcf"
        if not os.path.exists(f"{args.out_dir}/{chromo}"):
            os.mkdir(f"{args.out_dir}/{chromo}")
        with open(filenm, 'w') as new_vcf:
            for line in header:
                new_vcf.write(line)
            for line in all_lines:
                if line.split()[0] == chromo:
                    new_vcf.write(line)

if __name__ == '__main__':
    main()
