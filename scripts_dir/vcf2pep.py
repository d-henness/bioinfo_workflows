import argparse
import pyensembl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--release', type = int, help = 'Ensembl release')
    args = parser.parse_args()

    data = pyensembl.EnsemblRelease(args.release)

    genes = data.gene_ids_at_locus(1, 10361726)

if __name__ == '__main__':
    main()
