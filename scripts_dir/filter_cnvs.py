import sys
filenm1 = sys.argv[1]
filenm2 = sys.argv[2]
if (filenm1 == '') or (filenm2 == ''):
    sys.exit('Missing arguments')
with open(filenm1, 'r') as data:
    all_lines = data.readlines()
    with open(filenm2, 'w') as newfile:
        for i, line in enumerate(all_lines):
            if i == 0:
                newfile.write(line)
            elif line.split()[1] != line.split()[2]:
                newfile.write(line.strip('chr'))
