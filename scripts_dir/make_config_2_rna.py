import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser("Make config file for snakemake workflows. Directories containing the string 'tmp' anywhere in their name will automatically be excluded. The user will need to manually match the directories as tumour: normal pairs")
    parser.add_argument('-p', '--path_to_data', nargs = '+', required = True, help = """path to directory containing data""")
    parser.add_argument('dir_string', nargs = '+', help = """One or more strings that match the directories containing .fq.gz files. Example: \'python3 /path/to/make_config_2.py RG MF\' will match the directories RG7_1, RG7_PBMC, and MF5_PBMC, but not the directory RG7_1_tmp""")
    args = parser.parse_args()

    libs = []
    merge_libs = []

    for path in args.path_to_data:
        for i in sorted(os.listdir(path)):
            for dir_name_short in args.dir_string:
                if os.path.isdir(f"{path}/{i}") and (dir_name_short in i) and not ('tmp' in i):
                    reads_list = []
                    for seq in sorted(os.listdir(f"{path}/{i}")):
                        if (seq.endswith('.fastq') or seq.endswith('.fq.gz')) and ('cut_u' not in seq):
                            reads_list.append(f"{path}/{i}/{seq}")
                    if len(reads_list) > 0:
                        merge_lib_string = f"{i}:  ["
                        for j in range(int(len(reads_list)/2)):
                            lib_string = "{:}_{:}:  ['{:}', '{:}']".format(i, j, reads_list[(j * 2)], reads_list[(j * 2) + 1])
                            libs.append(lib_string)
                            merge_lib_string = "{:}{:}_{:}, ".format(merge_lib_string, i, j)
                        merge_lib_string = "{:}]".format(merge_lib_string[:-2])
                        merge_libs.append(merge_lib_string)

    print("rna_libraries:")
    for i in range(len(libs)):
        print("  {:}".format(libs[i]))
    print('rna_merge_libs:')
    for i in range(len(merge_libs)):
        print("  {:}".format(merge_libs[i]))

if __name__ == '__main__':
    main()
