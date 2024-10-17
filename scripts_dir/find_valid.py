import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenms', nargs = '+')
    parser.add_argument('--max_stragglers', type = int, default = 30) # largest amount of cancer stragglers to be considered dormant
    parser.add_argument('--min_cells_for_early_cuttoff', type = int, default = 40000) # smallest amount of cancer cells to be considered an early loss
    parser.add_argument('--max_timesteps', type = int, default = 14391) # if a file has fewer than this many timesteps it exited early
    parser.add_argument('--min_timesteps_for_domancy', type = int, default = 30) # if cancer population goes lower than max_stragglers for this many timesteps the tumor is dormant
    args = parser.parse_args()

    valid_files = []
    invalid_files = []

    for filenm in args.filenms:
        cancer_cells_per_timestep = []
        immune_cells_per_timestep = []
        with open(filenm, 'r') as data:
            for line in data:
                if line.startswith("Cancer cells"):
                    cancer_cells_per_timestep.append(int(line.split('\t')[1]))
                elif line.startswith("Immune cells"):
                    immune_cells_per_timestep.append(int(line.split('\t')[1]))
            if (len(cancer_cells_per_timestep) == 0):
                invalid_files.append(filenm)
            elif (len(cancer_cells_per_timestep) < args.max_timesteps) and (cancer_cells_per_timestep[-1] > args.max_stragglers) and (cancer_cells_per_timestep[-1] < args.min_cells_for_early_cuttoff):
                invalid_files.append(filenm)
            else:
                valid_files.append(filenm)
                print(len(cancer_cells_per_timestep))

#    with open("invalid_files.txt", 'w') as data:
#        for filenm in invalid_files:
#            data.write(f'{filenm}\n')
#
#    # get the estimated time of dormancy for vaild files
#    for i, filenm in enumerate(valid_files):
#        with open(filenm, 'r') as data:
#            cancer_cells_per_timestep = []
#            immune_cells_per_timestep = []
#            timesteps_lower_than_max_stragglers = 0
#            time_of_dormancy = -1
#            tumor_went_dormant = False
#            for line in data:
#                if line.startswith("Cancer cells"):
#                    ncancer_cells = int(line.split('\t')[1])
#                    if ncancer_cells < args.max_stragglers:
#                        timesteps_lower_than_max_stragglers += 1
#                    else:
#                        timesteps_lower_than_max_stragglers = 0
#                        tumor_went_dormant = False
#                        time_of_dormancy = -1
#                    cancer_cells_per_timestep.append(int(line.split('\t')[1]))
#                elif line.startswith("Immune cells"):
#                    immune_cells_per_timestep.append(int(line.split('\t')[1]))
#                if (timesteps_lower_than_max_stragglers >= args.min_timesteps_for_domancy) and not tumor_went_dormant:
#                    tumor_went_dormant = True
#                    time_of_dormancy = len(cancer_cells_per_timestep)
#        starting_immune_cells = filenm.split('/')[1]
#        mean_oncoprotien = filenm.split('/')[2]
#        if i == 0:
#            with open("valid_results.tsv", 'w') as data:
#                data.write(f'filename\tmean_oncoprotien\tstarting_immune_cells\ttimesteps\tdorment\ttime_of_dormancy\tcancer_cells_at_end\timmune_cells_at_end\n')
#                data.write(f'{filenm}\t{mean_oncoprotien}\t{starting_immune_cells}\t{len(cancer_cells_per_timestep)}\t{tumor_went_dormant}\t{time_of_dormancy}\t{cancer_cells_per_timestep[-1]}\t{immune_cells_per_timestep[-1]}\n')
#        else:
#            with open("valid_results.tsv", "a") as data:
#                data.write(f'{filenm}\t{mean_oncoprotien}\t{starting_immune_cells}\t{len(cancer_cells_per_timestep)}\t{tumor_went_dormant}\t{time_of_dormancy}\t{cancer_cells_per_timestep[-1]}\t{immune_cells_per_timestep[-1]}\n')
#





if __name__ == "__main__":
    main()
