import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vaf', type = float, default = 0.05)
parser.add_argument('filenms', nargs = '+')
args = parser.parse_args()

for filenm in args.filenms:
    muts_above_vaf = 0
    with open(filenm) as data:
        for line in data:
            split_line = line.strip().split('\t')
            ref = split_line[3]
            alt = split_line[4]
            format_helper_field = split_line[-3].split(':')
            tumor_format_field = split_line[-1].split(':')

            #somatic mutations
            if len(ref) == 1 and len(alt) == 1:
                ref_index = format_helper_field.index(f'{ref}U')
                alt_index = format_helper_field.index(f'{alt}U')

                ref_counts = int(tumor_format_field[ref_index].split(',')[0])
                alt_counts = int(tumor_format_field[alt_index].split(',')[0])
            else:
                ref_index = format_helper_field.index('TAR')
                alt_index = format_helper_field.index('TIR')

                ref_counts = int(tumor_format_field[ref_index].split(',')[0])
                alt_counts = int(tumor_format_field[alt_index].split(',')[0])

            if ((ref_counts + alt_counts) > 0) and (alt_counts / (ref_counts + alt_counts)) >= args.vaf:
                muts_above_vaf += 1
    print(f'{filenm}\t{muts_above_vaf}')




