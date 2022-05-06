import argparse
import sys
import os

con_types = [
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "downstream_gene_variant",
        "non_coding_transcript_exon_variant",
        "frameshift_variant",
        "inframe_deletion",
        "stop_gained",
        "inframe_insertion",
        "upstream_gene_variant",
        "3_prime_UTR_variant",
        "protein_altering_variant",
        "intergenic_variant",
        "mature_miRNA_variant",
        "coding_sequence_variant",
        "non_coding_transcript_variant",
        "splice_donor_variant",
        "NMD_transcript_variant",
        "intron_variant",
        "stop_retained_variant",
        "5_prime_UTR_variant",
        "missense_variant",
        "stop_lost",
        "start_retained_variant",
        "synonymous_variant",
        "start_lost",
        "splice_acceptor_variant"
]

skip_con_types = [
        "non_coding_transcript_exon_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_variant",
        "intron_variant",
        "5_prime_UTR_variant",
        "synonymous_variant",
]

parser = argparse.ArgumentParser()
parser.add_argument('filenms', nargs = '+')
args = parser.parse_args()

print("sample\tSNP\tindel")
for filenm in args.filenms:
    with open(filenm, 'r') as data:
        snp = 0
        indel = 0
        for i, line in enumerate(data):
            skip = False
            split_line = line.strip().split('\t')
            muts = split_line[7].split('&')
            for mut in muts:
                if mut in skip_con_types:
                    skip = True
                    break
            if not skip:
                ref = split_line[2].split(',')
                alts = split_line[3].split(',')
                for alt in alts:
                    if len(alt) != len(ref):
                        indel += 1
                    elif alt != ref:
                        snp += 1
                    else:
                        sys.exit(f"unparsable mutation in {filenm} at line {i}: {ref}, {alt}")
        print(f"{os.path.basename(filenm)[:-len('_VEP_parsed.tsv')]}\t{snp}\t{indel}")
