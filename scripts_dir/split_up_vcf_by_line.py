import argparse


def find_next_split(test_site, chro, all_sites):
    low = 0
    high = len(all_sites)
    while(low <= high):
        mid = int((low + high) / 2)
        test_num = int(all_sites[mid])
        if (mid == high):
            break
        if (test_num > test_site):
            high = mid
        elif (test_num <= test_site):
            low = mid + 1
    return mid, test_num


parser = argparse.ArgumentParser(description = 'Split up vcf file into chunks. Chunks will be as close together in size as possible.')
parser.add_argument('vcf_file', type = str, nargs = 1, help = 'vcf file to split up')
parser.add_argument('prefix', type = str, nargs = 1, help = 'file prefix for each new vcf file')
parser.add_argument('valid_sites_dir', type = str, nargs = 1, help = 'directory containing valid sites files')
parser.add_argument('nlines', type = int, nargs = 1, help = 'minimum number of lines per file')
parser.add_argument('--dry_run', action = 'store_true', help = 'only say how many files get created')
args = parser.parse_args()

with open(args.vcf_file[0], 'r') as data:
    all_file = data.readlines()
    lines_of_header = 0
    for line in all_file:
        lines_of_header += 1
        if line[0:len('#CHROM')] == '#CHROM':
            break
    chrs = {}
    for line in all_file[lines_of_header:]:
        try:
            chrs[line.split()[0]].append(line)
        except:
            chrs[line.split()[0]] = []
            chrs[line.split()[0]].append(line)

files_writen = 0
for key in chrs:
    marker = 0
    with open(args.valid_sites_dir[0] + f"/valid_sites_{key}.txt", "r") as data:
        all_sites = data.readlines()
        while (True):
            start = marker
            if (start + args.nlines[0] - 1 < len(chrs[key])):
                test_site = int(chrs[key][marker + args.nlines[0] - 1].split()[1])
            else:
                with open(args.prefix[0] + f"_split_{files_writen}.vcf", "w") as new_file:
                    new_file.writelines(all_file[0:lines_of_header])
                    new_file.writelines(chrs[key][start:])
                    files_writen += 1
                break
            next_split = find_next_split(test_site, key, all_sites)
            marker += args.nlines[0]

            while((marker < len(chrs[key])) and (int(chrs[key][marker].split()[1]) < next_split[1])):
                marker += 1
            if marker < len(chrs[key]):
                print(next_split[1], marker, chrs[key][marker - 1].split()[1], chrs[key][marker].split()[1])
                with open(args.prefix[0] + f"_split_{files_writen}.vcf", "w") as new_file:
                    new_file.writelines(all_file[0:lines_of_header])
                    new_file.writelines(chrs[key][start:marker])
                    files_writen += 1
            else:
                print(next_split[1], marker, chrs[key][marker - 1].split()[1])
                with open(args.prefix[0] + f"_split_{files_writen}.vcf", "w") as new_file:
                    new_file.writelines(all_file[0:lines_of_header])
                    new_file.writelines(chrs[key][start:])
                    files_writen += 1
                break



#marker = 0
#if not args.dry_run:
#    while (marker + lines_of_header < len(all_file)):
#        if (marker + args.nlines[0] - 1 < len(chr1)):
#            test_site = int(chr1[marker + args.nlines[0] - 1].split()[1])
#        else:
#            break
#
#        chro = chr1[marker + args.nlines[0] - 1].split()[0]
#        next_split = find_next_split(test_site, chro, args.valid_sites[0])
#        marker += args.nlines[0]
#
#        if (marker < len(chr1)):
#            while(int(chr1[marker].split()[1]) < next_split[1]):
#                marker += 1
#            print(next_split[1], marker, chr1[marker - 1].split()[1], chr1[marker].split()[1])
