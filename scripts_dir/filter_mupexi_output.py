import argparse
import operator

def parse_line(line, rna):
    parsed_line = {}
    line_split = line.strip().split('\t')
    # correct as of mupexi version 1.2.0
    mupexi_cols = ["HLA_allele", "Norm_peptide", "Norm_MHCrank_EL", "Mut_peptide", "Mut_MHCrank_EL", "Gene_ID", "Transcript_ID", "Amino_Acid_Change", "Allele_Frequency", "Mismatches", "peptide_position", "Chr", "Genomic_Position", "Protein_position", "Mutation_Consequence", "Gene_Symbol", "Cancer_Driver_Gene", "Proteome_Peptide_Match", "Expression_Level", "Mutant_affinity_score", "Normal_affinity_score", "Expression_score", "priority_Score"]
    for i, col in enumerate(mupexi_cols):
        if (i == 2) or (i == 4) or (i == 8):
                parsed_line[col] = float(line_split[i])
        elif rna and (i == 18):
                parsed_line[col] = float(line_split[i])
        else:
            parsed_line[col] = line_split[i]
    return parsed_line


parser = argparse.ArgumentParser(description = 'Filter MuPeXI output files for most relavent neoantigens')
parser.add_argument('mupexi_file', type = str, help = 'MuPeXI file to filter')
parser.add_argument('--rna', action = 'store_true', help = 'Toggle RNA expression filter')
args = parser.parse_args()

header = []
antigens = {}
whole_line = {}
with open(args.mupexi_file, 'r') as data:
    for i, line in enumerate(data):
        if i > 5:
            try:
                parsed_line = parse_line(line, args.rna)
            except:
                print(args.mupexi_file, line)
                quit()
            # Filter for neoantigens with Mut_MHC_rank_EL =< 0.5 AND Mut_peptide with exactly 9 letters AND Allele_frequency > 0.25
            if (parsed_line['Mut_MHCrank_EL'] <= 0.5) and (len(parsed_line['Mut_peptide']) == 9) and (parsed_line['Allele_Frequency'] > 0.25):
                if args.rna:
                    try:
                        if (antigens[parsed_line['Mut_peptide']] > parsed_line['Mut_MHCrank_EL']) and (parsed_line['Expression_Level'] > 0.1):
                            antigens[parsed_line['Mut_peptide']] = parsed_line['Mut_MHCrank_EL']# take the mut_peptide colunm as a key and figure out if the new prediction for that peptide has a lower binding strength then the current best prediction
                            whole_line[parsed_line['Mut_peptide']] = line
                    except:
                        if parsed_line['Expression_Level'] > 0.1:
                            antigens[parsed_line['Mut_peptide']] = parsed_line['Mut_MHCrank_EL']
                            whole_line[parsed_line['Mut_peptide']] = line
                else:
                    try:
                        if antigens[parsed_line['Mut_peptide']] > parsed_line['Mut_MHCrank_EL']:
                            antigens[parsed_line['Mut_peptide']] = parsed_line['Mut_MHCrank_EL']# take the mut_peptide colunm as a key and figure out if the new prediction for that peptide has a lower binding strength then the current best prediction
                            whole_line[parsed_line['Mut_peptide']] = line
                    except:
                        antigens[parsed_line['Mut_peptide']] = parsed_line['Mut_MHCrank_EL']
                        whole_line[parsed_line['Mut_peptide']] = line
        else:
            header.append(line.strip())

# print output to screen
for line in header:
    print(line)
bind_strength = None
for i, key_value in enumerate(sorted(antigens.items(), key = lambda kv: kv[1])):
#    print(f"{whole_line[key]}".strip())
    if (key_value[1] < 0.05) and (bind_strength != 'strong'):
        print("Strong binders")
        bind_strength = 'strong'
    elif (0.05 <= key_value[1] <= 0.15) and (bind_strength != 'inter'):
        print("Intermediate binders")
        bind_strength = 'inter'
    elif (key_value[1] > 0.15) and (bind_strength != 'weak'):
        print("Weak binders")
        bind_strength = 'weak'
    print(whole_line[key_value[0]])
