# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 16:04:45 2016

@author: dgreen
"""

# A crude plot of the 3D topography on Ascension?

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pandas as pd

# Read in the xyz file
dirpath = '/Users/dgreen/Documents/Work/4codor/topo_input_files/ascension/'
filename = 'ascension.xyz'
filename2 = 'temp_L1_L2_L3_L4_windfarm.dat'
df = pd.read_csv(dirpath+filename,skiprows=None,names=['x','y','z'],sep=r"\s*")
dfstats = pd.read_csv(dirpath+filename2,skiprows=None,names=['x','y'],sep=r"\s*")


para = {'axes.labelsize': 16, 'text.fontsize': 16, 'legend.fontsize': 14, 'xtick.labelsize': 14,'ytick.labelsize': 14, 'figure.subplot.left': 0.12, 'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11, 'figure.subplot.top': 0.97}
plt.rcParams.update(para)


fig = plt.figure(figsize=(5.0,5.5))
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_axes([0.25, 0.25, 0.75, 0.75])
ax.tricontourf(df['x'],df['y'],df['z'],40)

#Station
ax.plot(dfstats.iloc[0]['x'],dfstats.iloc[0]['y'],'ko')
ax.plot(dfstats.iloc[1]['x'],dfstats.iloc[1]['y'],'ko')
ax.plot(dfstats.iloc[2]['x'],dfstats.iloc[2]['y'],'ko')
ax.plot(dfstats.iloc[3]['x'],dfstats.iloc[3]['y'],'ko')
# Windfarm
ax.plot(dfstats.iloc[4]['x'],dfstats.iloc[4]['y'],'r*',markersize=10)
ax.set_aspect('equal', 'datalim')
ax.set_xlim([0,4000])



plt.show()
fig.savefig('ascension_xyz_example_crude.png',bbox_inches='tight')
