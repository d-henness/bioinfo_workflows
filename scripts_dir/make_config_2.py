import os
import sys

def main():
    try:
        dir_string = sys.argv[1]
    except:
        sys.exit("No dir string given")
    if dir_string == '':
        sys.exit("No dir string given")

    libs = []
    merge_libs = []

    for i in sorted(os.listdir('.')):
        if os.path.isdir(i) and (dir_string in i) and not ('tmp' in i):
            reads_list = []
            merge_lib_string = "{:}:  [".format(i)
            for seq in sorted(os.listdir(i)):
                if (('.fastq' in seq) or ('.gz' in seq)) and ('cut_u' not in seq):
                    reads_list.append(os.getcwd() + '/' + i + '/' + seq)
            for j in range(int(len(reads_list)/2)):
                lib_string = "{:}_{:}:  ['{:}', '{:}']".format(i, j, reads_list[(j * 2)], reads_list[(j * 2) + 1])
                libs.append(lib_string)
                merge_lib_string = "{:}{:}_{:}, ".format(merge_lib_string, i, j)
            merge_lib_string = "{:}]".format(merge_lib_string[:-2])
            merge_libs.append(merge_lib_string)

    print("libraries:")
    for i in range(len(libs)):
        print("  {:}".format(libs[i]))
    print('merge_libs:')
    for i in range(len(merge_libs)):
        print("  {:}".format(merge_libs[i]))

if __name__ == '__main__':
    main()
