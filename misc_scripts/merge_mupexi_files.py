import os

all_data = []
base_string = ''
# get data from all files with a mupexi extension
for directory in sorted(os.listdir('.')):
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
        else:
            print('Did not find mupexi file in {:}'.format(directory))

# write data to one file
with open(base_string[:len(base_string) - len('mupexi')] + 'merged.mupexi', 'w') as data:
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
