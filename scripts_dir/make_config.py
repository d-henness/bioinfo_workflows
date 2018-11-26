import os
import sys

def main():
    try:
        dir_string = sys.argv[1]
    except:
        sys.exit("No dir string given")
    if dir_string == '':
        sys.exit("No dir string given")
    reads_list = []
    dir_list = []
    for i in sorted(os.listdir('.')):
        if os.path.isdir(i) and (dir_string in i) and not ('tmp' in i):
            dir_list.append(i)
            for seq in sorted(os.listdir(i)):
                if (('.fastq' in seq) or ('.gz' in seq)) and ('cut_u' not in seq):
                    reads_list.append(os.getcwd() + '/' + i + '/' + seq)

    file_pairs = []
    for i in range(0, len(reads_list), 2):
        file_pairs.append((reads_list[i], reads_list[i + 1]))

    j = 0
    k = -1
    lib_name = dir_list[j]
    libs = []
    merge_libs = ["{:}: [".format(lib_name)]

    for i in range(len(file_pairs)):
        if lib_name in file_pairs[i][0]:
            k += 1
            merge_libs[j] = "{:}{:}_{:},".format(merge_libs[j], lib_name, k)
        else:
            merge_libs[j] = "{:}]".format(merge_libs[j][:-1]) # remove trailing comma
            k = 0
            j += 1
            lib_name = dir_list[j]
            merge_libs.append("{:}: [{:}_{:},".format(lib_name, lib_name, k))
        libs.append("{:}_{:}: ['{:}', '{:}']".format(lib_name, k, file_pairs[i][0], file_pairs[i][1]))
    merge_libs[-1] = merge_libs[-1][:-1] + ']'

    print('libraries:')
    for i in range(len(libs)):
        print("  {:}".format(libs[i]))
    print('merge_libs:')
    for i in range(len(merge_libs)):
        print("  {:}".format(merge_libs[i]))


if __name__ == '__main__':
    main()
