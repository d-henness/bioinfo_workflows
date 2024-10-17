import argparse
import numpy as np
import os
import matplotlib.pyplot as plt

def make_onco_hist(args, fig, ax):
    oncoprotien_all_timesteps = []
    oncoprotien_this_timestep = []
    copying_data = False

    for i, filenm in enumerate(args.filenms):
        print(filenm)
        timestep = 0
        with open(filenm, 'r') as data:
            for line in data:
                if line.startswith('Oncoprotien\tposition'):
                    copying_data = True
                elif line.startswith('Cancer cells:'):
                    if i == 0:
                        oncoprotien_all_timesteps.append(oncoprotien_this_timestep)
                    else:
                        if len(oncoprotien_all_timesteps) > timestep:
                            oncoprotien_all_timesteps[timestep] += oncoprotien_this_timestep
                        else:
                            oncoprotien_all_timesteps.append(oncoprotien_this_timestep)
                    oncoprotien_this_timestep = []
                    copying_data = False
                    timestep += 1
                elif copying_data:
                    # check if cell is stragller
                    split_line = line.strip().split('\t')
                    pos = split_line[1].split(' ')
                    stragller = False
                    for coord in pos:
                        if abs(float(coord)) > 500: # change if using a different area than 500 x 500
                            stragller = True
                            break
                    if not stragller:
                        oncoprotien_this_timestep.append(float(split_line[0]))

        # Flatten the data for the histogram and create corresponding timestep and index arrays
    x_data = []  # timestep
    y_data = []  # Oncoprotien

    for i, timestep in enumerate(oncoprotien_all_timesteps):
        for _, oncoprotien in enumerate(timestep):
            x_data.append(i)
            y_data.append(oncoprotien)

    x_data = np.array(x_data)
    y_data = np.array(y_data)


    x_bins = np.array([i * 10 for i in range(int(len(oncoprotien_all_timesteps)/10) + 1)])
    y_bins = np.array([i * 0.1 for i in range(21)])
    hist, xedges, yedges = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins])
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    dx = (xedges[1] - xedges[0]) * 0.9 * np.ones_like(zpos)
    dy = (yedges[1] - yedges[0]) * 0.9 * np.ones_like(zpos)
    dz = hist.ravel()
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
    ax.set_xlabel('Time Steps')
    ax.set_ylabel('Oncoprotien')
    ax.set_zlabel('Cancer Cell Counts')


    #get the get the start, min, and max after min for each set of oncoprotien bins
    stats_filenm = os.path.splitext(args.outfile)[0] + ".tsv"
    with open(stats_filenm, 'w') as data:
        data.write('onco\tstart\tmin\tmax_after_min\n')
        for i, onco_data in enumerate(hist.T):
            start = onco_data[0]
            min_o = min(onco_data)
            min_o_index = list(onco_data).index(min_o)
            max_after_min = max(onco_data[min_o_index:])
            data.write(f'{i * 0.1:0.2f}\t{start}\t{min_o}\t{max_after_min}\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenms', nargs = '+')
    parser.add_argument('-e', '--elevation', type = float, default = 30)
    parser.add_argument('-a', '--azimuth', type = float, default = 45)
    parser.add_argument('-o', '--outfile', default = 'out.pdf')
    args = parser.parse_args()

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.view_init(elev = args.elevation, azim = args.azimuth)

    make_onco_hist(args, fig, ax)
    fig.tight_layout()
    plt.savefig(args.outfile)

if __name__ == '__main__':
    main()
