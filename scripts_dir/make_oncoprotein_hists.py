import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenms', nargs = '+')
    args = parser.parse_args()

    oncoprotien_all_timesteps = []
    oncoprotien_this_timestep = []
    copying_data = False

    for i, filenm in enumerate(args.filenms):
        timestep = 0
        with open(filenm, 'r') as data:
            for line in data:
                if line.startswith('Oncoprotien\tposition'):
                    copying_data = True
                elif line.startswith('Cancer cells:'):
                    if i == 0:
                        oncoprotien_all_timesteps.append(oncoprotien_this_timestep)
                    else:
                        oncoprotien_all_timesteps[i] += oncoprotien_this_timestep
                    oncoprotien_this_timestep = []
                    copying_data = False
                    timestep += 1
                elif copying_data:
                    oncoprotien_this_timestep.append(float(line.split('\t')[0]))
        print(oncoprotien_all_timesteps)

        # Flatten the data for the histogram and create corresponding timestep and index arrays
    x_data = []  # timestep
    y_data = []  # Oncoprotien

    for i, timestep in enumerate(oncoprotien_all_timesteps):
        for j, oncoprotien in enumerate(timestep):
            x_data.append(i)
            y_data.append(oncoprotien)

    x_data = np.array(x_data)
    y_data = np.array(y_data)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    hist, xedges, yedges = np.histogram2d(x_data, y_data, bins=[int(len(oncoprotien_all_timesteps)/10), 20])
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    print(xpos)
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    dx = (xedges[1] - xedges[0]) * 0.9 * np.ones_like(zpos)
    dy = (yedges[1] - yedges[0]) * 0.9 * np.ones_like(zpos)
    dz = hist.ravel()
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')

    plt.savefig("test.pdf")

if __name__ == '__main__':
    main()
