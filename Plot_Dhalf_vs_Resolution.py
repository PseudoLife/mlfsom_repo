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
image_folder = '/home/peter/Desktop/MLFSOM/data_homo_mosOff_1.2A_10fr'

#df_peaks = myfunc.ReadGaussianPeakIntensities(image_folder)
df_peaks = myfunc.ReadPeakIntensities(image_folder)
df_shells = myfunc.ShellIntensities(df_peaks,20)

plt.close('all')
fig, ax = plt.subplots()
plt.ion()

# data
d_space = df_shells.d_max*0.5+df_shells.d_min*0.5
D_half = df_shells.D_half
plt.scatter(d_space,D_half,edgecolors='black',s=150)
# fitting line
x = np.linspace(1.2, 5.4, 100)
K = 0.36
y = K * x**2
ax.plot(x, y, linestyle = '--', color='tab:red', linewidth=2)


#ax.set_ylim(8,100)
#ax.set_yticks([10,100])
#ax.set_xlim(1,10)
#ax.set_xticks([1,10])
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

#fig.savefig('fig_DhalfvsRes',dpi=300)
plt.show()