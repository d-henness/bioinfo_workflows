import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenm')
    args = parser.parse_args()

    with open(args.filenm, 'r') as data:
        cols = {}
        for i, line in enumerate(data):
            split_line = line.strip().split('\t')
            if i == 0:
                for j, col in enumerate(split_line):
                    cols[col] = j
                print(f'CHROM\tPOS\tREF\tALT\tVAF\tSYMBOL\tGene\tConsequence\tSIFT\tPolyPhen\tCondel\tLoFtool\tBLOSUM62\tProtein_position\tAmino_acids\tCodons')
            else:
                # snps
                alt_counts = None
                ref_counts = None
                if len(split_line[cols['REF']]) == 1 and len(split_line[cols['ALT']]) == 1:
                    alt_counts = int(split_line[cols[split_line[cols['ALT']]]].split(',')[0])
                    ref_counts = int(split_line[cols[split_line[cols['REF']]]].split(',')[0])
                else:
                    alt_counts = int(split_line[cols['TIR']].split(',')[0])
                    ref_counts = int(split_line[cols['TAR']].split(',')[0])
                if (alt_counts + ref_counts != 0):
                    VAF = alt_counts / (alt_counts + ref_counts)
                    new_line = f'{split_line[cols["CHROM"]]}\t{split_line[cols["POS"]]}\t{split_line[cols["REF"]]}\t{split_line[cols["ALT"]]}\t{VAF}\t{split_line[cols["SYMBOL"]]}\t{split_line[cols["Gene"]]}\t{split_line[cols["Consequence"]]}\t{split_line[cols["SIFT"]]}\t{split_line[cols["PolyPhen"]]}\t{split_line[cols["Condel"]]}\t{split_line[cols["LoFtool"]]}\t{split_line[cols["BLOSUM62"]]}\t{split_line[cols["Protein_position"]]}\t{split_line[cols["Amino_acids"]]}\t{split_line[cols["Codons"]]}'
                    print(new_line)
                else:
                    print(f'problem line {line} skipping', file = sys.stderr)






if __name__ == '__main__':
    main()
