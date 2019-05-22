import argparse
import operator

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
            split_line = line.split()
            if float(split_line[4]) <= 0.5:
                if args.rna:
                    try:
                        if (antigens[split_line[3]] > float(split_line[4])) and (float(split_line[18]) > 0.1):
                            antigens[split_line[3]] = float(split_line[4])# take the mut_peptide colunm as a key and figure out if the new prediction for that peptide has a lower binding strength then the current best prediction
                            whole_line[split_line[3]] = line
                    except:
                        if float(split_line[18]) > 0.1:
                            antigens[split_line[3]] = float(split_line[4])
                            whole_line[split_line[3]] = line
                else:
                    try:
                        if antigens[split_line[3]] > float(split_line[4]):
                            antigens[split_line[3]] = float(split_line[4])# take the mut_peptide colunm as a key and figure out if the new prediction for that peptide has a lower binding strength then the current best prediction
                            whole_line[split_line[3]] = line
                    except:
                        antigens[split_line[3]] = float(split_line[4])
                        whole_line[split_line[3]] = line
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
