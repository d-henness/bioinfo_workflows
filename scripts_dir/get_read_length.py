import sys
from xopen import xopen


def main():
    read1 = sys.argv[1]

    with xopen(read1, 'r') as gatk_file:
        line = gatk_file.readline()
        line = gatk_file.readline()
        length = len(line)
        print(length)


if __name__ == '__main__':
    main()
