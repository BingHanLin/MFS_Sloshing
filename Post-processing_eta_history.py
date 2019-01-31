import os
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

plt.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath} \boldmath']
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rc('font', weight='bold')

Nt = 30
dt = 0.05
offest = 0.3
fig, ax = plt.subplots(figsize=(6, 4), dpi=300, facecolor='w', edgecolor='k')

#============================================================
# Read data from solution data
#============================================================
dirname = 'OutputData'
y1_list = []
y2_list = []
time_list = []

Time = 0

for j in range(len(os.listdir(dirname))):
    Time = (j+1)*dt
    filename = '%0.5f.dat' % Time
    print (filename)

    filename = './' + dirname + '/' + filename

    dayatype = np.dtype({'names':['Nbou_x', 'Nbou_y', 'Nsou_x', 'Nsou_y'], 
                         'formats':['f', 'f','f', 'f']})

    data = np.loadtxt(filename, dtype = dayatype, delimiter=',')

    y1_list.append( (data['Nbou_y'][0]-offest    ) )
    y2_list.append( (data['Nbou_y'][Nt-1]-offest ) )

    time_list.append( Time )

line2 = ax.plot(time_list, y2_list, linewidth = 1, color="black", label = 'Present result')

#============================================================
# Read data from example
#============================================================
filename = 'BFC_example.txt'

dayatype = np.dtype({'names':['time', 'Eta'], 'formats':['f', 'f']})

data_example = np.loadtxt(filename, dtype = dayatype, delimiter=None)

line3 = ax.plot(data_example['time'], data_example['Eta']/100,linestyle='-.', 
                linewidth = 1, color="red", label = 'BFC et al.')

#============================================================
# Figure
#============================================================
ax.set_xlim(0,  9)
ax.set_xlabel(r'$\textrm{Time (sec)}$', fontsize=15, labelpad=10)
ax.set_ylabel(r'$\textrm{Wave Height (m)}$', fontsize=15, labelpad=10)
plt.legend( prop={'size': 10}, ncol=1, edgecolor = 'k')
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig('output.png')
plt.show()