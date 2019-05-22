import os
import argparse

parser = argparse.ArgumentParser(description = 'merge mupexi data')
parser.add_argument('base_dir', type = str, help = 'directory to look for mupexi files')
parser.add_argument('prefix', type = str, help = 'prefix for output files')
args = parser.parse_args()

all_data = []
base_string = ''
missense = 0
insertions = 0
deletions = 0
frameshift = 0
# get data from all files with a mupexi extension
for directory in sorted(os.listdir(args.base_dir)):
    if os.path.isdir(directory):
        found = False
        mupexi_file = ''
        for filenm in sorted(os.listdir(directory)):
            if (filenm[len(filenm) - 6:len(filenm)] == 'mupexi') and (os.path.isfile(directory + '/' + filenm)):
                found = True
                mupexi_file = filenm
                base_string = filenm
                break
        if found:
            with open(directory + '/' + mupexi_file) as data:
                all_data.append(data.readlines())
            mupexi_log = mupexi_file[:-len(".mupexi")] + ".log"
            with open(directory + '/' + mupexi_log) as data:
                all_lines = data.readlines()
                for i, line in enumerate(all_lines):
                    if "missense" in line:
                        missense += int(line.split()[1])
                        insertions += int(all_lines[i+1].split()[0])
                        deletions += int(all_lines[i+2].split()[0])
                        frameshift += int(all_lines[i+3].split()[0])
                        break
        else:
            print('Did not find mupexi file in {:}'.format(directory))

# write data to one file
with open(f"{args.base_dir}/{args.prefix}_merged.mupexi", 'w') as data:
    for i, file_data in enumerate(all_data):
        if i == 0:
            data.write(''.join(file_data))
        else:
            # find the length of the header
            header_len = 0
            for line in file_data:
                if line[:len('HLA-')] == 'HLA-':
                    break
                else:
                    header_len += 1
            data.write(''.join(file_data[header_len:]))

with open(f"{args.base_dir}/{args.prefix}_merged.log", 'w') as data:
    data.write("missense {}\n".format(missense))
    data.write("insertions {}\n".format(insertions))
    data.write("deletions {}\n".format(deletions))
    data.write("frameshift {}\n".format(frameshift))
