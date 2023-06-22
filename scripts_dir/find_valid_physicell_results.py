import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenms', nargs = '+')
    parser.add_argument('--min_cells_for_early_cuttoff', type = int, default = 40000) # smallest amount of cancer cells to be considered an early loss
    parser.add_argument('--max_timesteps', type = int, default = 14391) # if a file has fewer than this many timesteps it exited early
    parser.add_argument('--invalid_files', '-i', nargs = '+')
    args = parser.parse_args()

    valid_files = []
    invalid_files = []

    with open('valid_files.tsv', 'w') as valid_file_data:
        valid_file_data.write(f'filename\tmean_oncoprotien_at_start\tstarting_immune_cells\ttimesteps\tdorment\ttime_of_dormancy\tcancer_cells_at_end\timmune_cells_at_end\n')
    with open('invalid_files.txt', 'w') as invalid_file_data:
        invalid_file_data.write('invalid_files\n')
        if args.invalid_files is not None:
            invalid_file_data.writelines(args.invalid_files)

    for filenm in args.filenms:
        if filenm not in invalid_files:
            print(filenm)
            timestep = 0
            cancer_cells_all_timesteps = []
            cancer_cells_this_timestep = 0
            counting_cells = False
            eradicated = False
            valid = False
            time_of_eradication = -1
            immune_cells_at_end = 0

            with open(filenm, 'r') as data:
                for line in data:
                    if line.startswith('Oncoprotien\tposition'):
                        counting_cells = True
                    elif line.startswith('Immune cells:'):
                        immune_cells_at_end = int(line.split('\t')[1]) # this will be writen over each time so at end of the loop it will be the last count
                    elif line.startswith('Cancer cells:'):
                        cancer_cells_all_timesteps.append(cancer_cells_this_timestep)
                        if (cancer_cells_this_timestep == 0) and not eradicated:
                            eradicated = True
                            time_of_eradication = timestep
                        cancer_cells_this_timestep = 0
                        counting_cells = False
                        timestep += 1
                    elif counting_cells:
                        # check if cell is stragller
                        split_line = line.strip().split('\t')
                        pos = split_line[1].split(' ')
                        stragller = False
                        if len(pos) == 3:
                            for coord in pos:
                                if abs(float(coord)) > 500: # change if using a different area than 500 x 500
                                    stragller = True
                                    break
                            if not stragller:
                                cancer_cells_this_timestep += 1

            print(timestep, filenm)
            if eradicated or (timestep >= args.max_timesteps):
                valid_files.append(filenm)
                valid = True
            elif len(cancer_cells_all_timesteps) == 0:
                invalid_files.append(filenm)
            elif cancer_cells_all_timesteps[-1] >= args.min_cells_for_early_cuttoff:
                valid_files.append(filenm)
                valid = True
            else:
                invalid_files.append(filenm)

            if valid:
                with open('valid_files.tsv', 'a') as valid_file_data:
                    mean_onco = filenm.split('/')[2]
                    starting_immune = filenm.split('/')[1]
                    valid_file_data.write(f'{filenm}\t{mean_onco}\t{starting_immune}\t{timestep - 1}\t{eradicated}\t{time_of_eradication}\t{cancer_cells_all_timesteps[-1]}\t{immune_cells_at_end}\n')
            else:
                with open('invalid_files.txt', 'a') as invalid_file_data:
                    invalid_file_data.write(f'{filenm}\n')



if __name__ == '__main__':
    main()
