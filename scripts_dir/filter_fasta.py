import argparse

def main():
    parser = argparse.ArgumentParser(description = 'filter fasta file for only the primary regions')
    parser.add_argument('fasta_file', type = str, help = 'fasta file to filter')
    args = parser.parse_args()

    chrs = [
        'chr1',
        'chr2',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr7',
        'chr8',
        'chr9',
        'chr10',
        'chr11',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        'chr17',
        'chr18',
        'chr19',
        'chr20',
        'chr21',
        'chr22',
        'chrX',
        'chrY'
    ]

    with open(args.fasta_file, 'r') as data:
        should_write = False
        for line in data:
            if (line[0] == '>'):
                if (line.split()[0][1:] in chrs):
                    should_write = True
                else:
                    should_write = False
            if should_write:
                print(line, end = '')

if __name__ == '__main__':
    main()
