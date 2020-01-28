import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser("Make config file for snakemake workflows. Directories containing the string 'tmp' anywhere in their name will automatically be excluded. The user will need to manually match the directories as tumour: normal pairs")
    parser.add_argument('dir_string', nargs = '+', help = """One or more strings that match the directories containing .fq.gz files. Example: \'python3 /path/to/make_config_2.py RG MF\' will match the directories RG7_1, RG7_PBMC, and MF5_PBMC, but not the directory RG7_1_tmp""")
    parser.add_argument('--bed_file', help = "bed_file for this run", required = True)
    args = parser.parse_args()

    libs = []
    merge_libs = []


    for dir_name_short in args.dir_string:
        path, basename = os.path.split(dir_name_short)
        print(path, basename)
        for i in sorted(os.listdir(path)):
            if (os.path.isdir(f"{path}/{i}") or os.path.islink(f"{path}/{i}")) and (i.startswith(basename)) and not ('tmp' in i):
                reads_list = []
                merge_lib_string = "{:}:  [".format(i)
                for seq in sorted(os.listdir(f"{path}/{i}")):
                    if (seq.endswith('.fastq') or seq.endswith('.fq.gz')) and ('cut_u' not in seq):
                        reads_list.append(os.getcwd() + '/' + path + '/' + i + '/' + seq)
                for j in range(int(len(reads_list)/2)):
                    lib_string = "{:}_{:}:  ['{:}', '{:}']".format(i, j, reads_list[(j * 2)], reads_list[(j * 2) + 1])
                    libs.append(lib_string)
                    merge_lib_string = "{:}{:}_{:}, ".format(merge_lib_string, i, j)
                merge_lib_string = "{:}]".format(merge_lib_string[:-2])
                merge_libs.append(merge_lib_string)

    print(f"alt_bed:  {args.bed_file}")

    print("libraries:")
    for i in range(len(libs)):
        print("  {:}".format(libs[i]))
    print('merge_libs:')
    for i in range(len(merge_libs)):
        print("  {:}".format(merge_libs[i]))

if __name__ == '__main__':
    main()
