import pandas as pd
import os, pickle
from os.path import join
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from matplotlib.ticker import AutoMinorLocator

from toPrecision import toPrecision
import FunctionsMLFSOMPeaks as myfunc

#image_folder = '/Users/atakisi/Desktop/MLFSOM/data_gaussian_N6_1.5A'
#image_folder = '/Users/atakisi/Desktop/MLFSOM/data_homo_mosOff_finer_Newbie'
#image_folder = '/home/peter/Desktop/MLFSOM/data_homo_mosOff_1.2A_10fr'

#df_peaks = myfunc.ReadGaussianPeakIntensities(image_folder)
#df_peaks = myfunc.ReadPeakIntensities(image_folder)
#df_shells = myfunc.ShellIntensities(df_peaks,20)
df_shells = pickle.load(open('/home/thorne20/Desktop/MLFSOM/data_stills_2.0A_50fr_1deg/xds_shells.pickle'))
K = 1.40
weights = (0.75,0.25) # weights of d_min, d_max in d_space = d_min*w_min + d_max*w_max

plt.close('all')
fig, ax = plt.subplots()
plt.ion()

# data
d_space = df_shells.d_min*weights[0] + df_shells.d_max*weights[1]
D_half = df_shells.D_half
plt.scatter(d_space,D_half,edgecolors='black',s=150)
# fitting line
x_min = d_space.iloc[-1] * 0.90
x_max = d_space.iloc[0] * 1.1
x = np.linspace(x_min, x_max, 100)
y = K * x**2
ax.plot(x, y, linestyle = '--', color='tab:red', linewidth=2)

ax.annotate(' $D_{1/2} \propto d^{2}$', xy=(0.74, 0.61), xycoords='axes fraction',\
	fontsize='xx-large', horizontalalignment='center', verticalalignment='center', color='tab:red')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Resolution ($\mathrm{\\AA}$)',fontsize='x-large')
ax.set_ylabel('Half-dose (arbitrary)',fontsize='x-large')
ax.tick_params(axis='both', which='both', length=0, labelsize='x-large',pad=10)
ax.grid(which='major',linewidth=0.4)
ax.grid(which='minor',linewidth=0.2)
ax.set_facecolor('0.95')
for axis in ['top','bottom','left','right']: ax.spines[axis].set_visible(False)
plt.tight_layout()

#fig.savefig('fig_DhalfvsRes',dpi=200)
plt.show()