import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenm')
    parser.add_argument('--min_cells_for_early_cuttoff', type = int, default = 40000) # smallest amount of cancer cells to be considered an early loss
    parser.add_argument('--max_timesteps', type = int, default = 14391) # if a file has fewer than this many timesteps it exited early
    parser.add_argument('--invalid_files', '-i', nargs = '+')
    args = parser.parse_args()

    valid_files = []
    invalid_files = []

    valid_results_file = f'{os.path.split(args.filenm)[0]}/valid_results.tsv'
    with open(valid_results_file, 'w') as valid_file_data:
        valid_file_data.write(f'filename\tnum_starting_cancer\tstarting_immune_cells\tmean_oncoprotien_at_start\trandom_seed\ttimestep\toncoprotien\n')

        print(args.filenm)
        timestep = 0
        cancer_cells_all_timesteps = []
        cancer_cells_this_timestep = 0
        counting_cells = False
        eradicated = False
        valid = False

        split_filenm = args.filenm.split('/')
        starting_cancer = split_filenm[0]
        starting_immune = split_filenm[1]
        starting_onco = split_filenm[2]
        random_seed = split_filenm[3]

        with open(args.filenm, 'r') as data:
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
                            valid_file_data.write(f'{args.filenm}\t{starting_cancer}\t{starting_immune}\t{starting_onco}\t{random_seed}\t{timestep}\t{split_line[0].strip()}\n')
                            cancer_cells_this_timestep += 1

    print(timestep, args.filenm)
    if eradicated or (timestep >= args.max_timesteps):
        valid_files.append(args.filenm)
        valid = True
    elif len(cancer_cells_all_timesteps) == 0:
        invalid_files.append(args.filenm)
    elif cancer_cells_all_timesteps[-1] >= args.min_cells_for_early_cuttoff:
        valid_files.append(args.filenm)
        valid = True
    else:
        invalid_files.append(args.filenm)

    if not valid:
        with open(valid_results_file, 'w') as invalid_file_data:
            invalid_file_data.write(f'invalid file\n')



if __name__ == '__main__':
    main()
