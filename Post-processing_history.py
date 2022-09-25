import matplotlib.animation as animation
from matplotlib import pyplot as plt
import os
from matplotlib import style
from matplotlib import animation
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl
from matplotlib.animation import FuncAnimation
# plt.rc('text', usetex=True)
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath} \boldmath']
# mpl.rcParams['axes.linewidth'] = 1.5
# mpl.rc('font', weight='bold')


def plot_wave_profile():
    Nt = 30
    dt = 0.05
    offest = 0.3
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300,
                           facecolor='w', edgecolor='k')

    # ============================================================
    # Read data from solution data
    # ============================================================
    dirname = 'OutputData'
    y1_list = []
    y2_list = []
    time_list = []

    time = 0

    for j in range(len(os.listdir(dirname))):
        time = (j+1)*dt
        filename = '%0.5f.dat' % time
        print(filename)

        filename = './' + dirname + '/' + filename

        dayatype = np.dtype({'names': ['Nbou_x', 'Nbou_y', 'Nsou_x', 'Nsou_y'],
                            'formats': ['f', 'f', 'f', 'f']})

        data = np.loadtxt(filename, dtype=dayatype, delimiter=',')

        y1_list.append((data['Nbou_y'][0]-offest))
        y2_list.append((data['Nbou_y'][Nt-1]-offest))

        time_list.append(time)

    line2 = ax.plot(time_list, y2_list, linewidth=1,
                    color="black", label='Present result')

    # ============================================================
    # Read data from example
    # ============================================================
    filename = 'BFC_example.txt'

    dayatype = np.dtype({'names': ['time', 'Eta'], 'formats': ['f', 'f']})

    data_example = np.loadtxt(filename, dtype=dayatype, delimiter=None)

    line3 = ax.plot(data_example['time'], data_example['Eta']/100, linestyle='-.',
                    linewidth=1, color="red", label='BFC et al.')

    # ============================================================
    # Figure
    # ============================================================
    ax.set_xlim(0,  9)
    ax.set_xlabel(f'Time (sec)', fontsize=15, labelpad=10)
    ax.set_ylabel(f'Wave Height (m)', fontsize=15, labelpad=10)
    plt.legend(prop={'size': 10}, ncol=1, edgecolor='k')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig('wave_profile.png')
    plt.show()


def plot_wave_animation():
    Nt = 30
    dt = 0.05

    dirname = 'OutputData'

    fig, ax = plt.subplots(figsize=(6, 4), dpi=300,
                           facecolor='w', edgecolor='k')

    def update_animate(i):
        fig.clear()

        ax = fig.add_subplot(aspect='equal')

        time = (i+1)*dt
        filename = '%0.5f.dat' % time
        print(filename)

        filename = './' + dirname + '/' + filename

        dayatype = np.dtype({'names': ['Nbou_x', 'Nbou_y', 'Nsou_x', 'Nsou_y'],
                            'formats': ['f', 'f', 'f', 'f']})

        data = np.loadtxt(filename, dtype=dayatype, delimiter=',')

        surface_list_x = data['Nbou_x'][0:Nt]
        surface_list_y = data['Nbou_y'][0:Nt]

        bou_list_x = data['Nbou_x'][Nt:-1]
        bou_list_y = data['Nbou_y'][Nt:-1]

        ax.plot([-0.5, 0.5], [0.3, 0.3], linestyle='--', linewidth=0.5,
                color='black')

        ax.scatter(surface_list_x, surface_list_y, linewidth=1,
                   color='blue', s=2.0, label='Surface Nodes')

        ax.scatter(bou_list_x, bou_list_y, linewidth=1,
                   color='red', s=2.0, label='Wall Nodes')

        ax.set_xlim(-0.5,  0.5)
        ax.set_ylim(-0.4,  0.6)
        ax.set_xlabel(f'X (m)', fontsize=15, labelpad=10)
        ax.set_ylabel(f'Y (m)', fontsize=15, labelpad=10)

        plt.legend(prop={'size': 10}, ncol=1, edgecolor='k')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.grid(color='gray', linestyle='-', linewidth=0.25)

    fps = 24
    n_frames = len(os.listdir(dirname))
    print(n_frames, ".....")
    ani = FuncAnimation(fig, update_animate,
                        frames=n_frames, interval=1000/fps)

    ani.save('wave_anim.gif', writer='pillow')


if __name__ == "__main__":
    plot_wave_animation()
