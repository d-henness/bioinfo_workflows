import argparse
import sys

parser = argparse.ArgumentParser(description = "Parse the output from VEP vcf files into a more tractable format. All vcf files must have the same VEP fields present")
parser.add_argument("-f", "--fields", nargs = '+', default = ["all"], help = "the fields wanted in the final output")
parser.add_argument("-v", "--vcf", nargs = '+', help = "vcf files to run on", required = True)
args = parser.parse_args()

for i, filenm in enumerate(args.vcf):
    with open (filenm, 'r') as data:
        for line in data:
            if (i == 0) and (line[0:len("##INFO=<ID=CSQ")] == "##INFO=<ID=CSQ"):
                vep_fields_global = line.strip().split()[-1][:-len("\">")].split("|")
                if (args.fields[0] != "all"):
                    vep_field_map = {}
                    for field in args.fields:
                        if (field not in vep_fields_global):
                            sys.exit(f"{field} not in vep output\nvalid fields are {vep_fields_global}")
                        else:
                            vep_field_map[field] = vep_fields_global.index(field)
                    print(f"chr,pos,ref,alt,{','.join([field for field in args.fields])}")
                else:
                    print(f"chr,pos,ref,alt,{','.join([field for field in vep_fields_global])}")

            elif (line[0:len("chr")] == "chr"):
                split_line = line.strip().split()
                vep_fields_line = split_line[7].split(";")[-1].split("|")
                if (args.fields[0] != "all"):
                    print(f"{split_line[0]},{split_line[1]},{split_line[3]},{split_line[4]},{','.join([vep_fields_line[vep_field_map[field]] for field in args.fields])}")
                else:
                    print(f"{split_line[0]},{split_line[1]},{split_line[3]},{split_line[4]},{','.join([field for field in vep_fields_line])}")
