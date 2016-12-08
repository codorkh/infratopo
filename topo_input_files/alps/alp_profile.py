# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 14:36:32 2016

@author: dgreen
"""

# alp_profile.py

# PLotting the alpine profile, chosen to give a relatively 'up and over' profile
# from the coordinates 44.27N 10.60E

# Load in the data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_profile(filename):
        data = pd.io.parsers.read_csv(filename, sep=r'\s*',names=['lon','lat','dist','alt'])
        return data
  
def read_2D_profile(filename):
        data = pd.io.parsers.read_csv(filename, sep=r'\s*',names=['dist','alt'])
        return data
      
para = {'axes.labelsize': 18, 'text.fontsize': 18, 'legend.fontsize': 13, 'xtick.labelsize': 16,'ytick.labelsize': 16, 'figure.subplot.left': 0.12, 'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11, 'figure.subplot.top': 0.97}
plt.rcParams.update(para)

dirpath = '/Users/dgreen/Documents/Work/4codor/infratopo/topo_input_files/alps/' 

profiledata = read_profile(dirpath+'alp.xydz')

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_axes([0.15,0.15,0.8,0.75])
ax1.plot(profiledata['dist'],profiledata['alt'],'k-')
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Elevation (m)')

fig.savefig(dirpath+'alpine_profile.png',bbox_inches='tight')

rzfile = 'alp_2d.dat'

fL1 = open(dirpath+rzfile,'w')
for x in range(len(profiledata['dist'])):
        fL1.write('{:8.1f} {:7.3f}\n'.format(profiledata.iloc[x]['dist']*1000.,profiledata.iloc[x]['alt']))
fL1.close()



fig = plt.figure(figsize=(10,5))
ax1 = fig.add_axes([0.15,0.15,0.8,0.75])
ax1.plot(profiledata['dist'],profiledata['alt'],'k-',label='Alps')
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Elevation (m)')

profilegaussdata = read_2D_profile('/Users/dgreen/Documents/Work/4codor/infratopo/topo_input_files/synthetics/gauss_3000m_hill_long.dat')

ax1.plot(profilegaussdata['dist']/1000.,profilegaussdata['alt'],'r-',label='Gaussian Synthetic')
ax1.legend(loc=1)
fig.savefig(dirpath+'alpine_profile_compare.png',bbox_inches='tight')