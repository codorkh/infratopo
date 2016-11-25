# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:42:12 2016

@author: dgreen
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Generating the input atmosphere files for Codor's PE code.
# Codor tells me he only needs the velocity profiles with altitude
# I will rip these out of the NCPA toy model for the test cases

#----------------------------------------------------------------
# Input data
dirstr = '/Users/dgreen/Documents/Work/4codor/atmos_input_files/'
filename = 'NCPA_canonical_profile_zuvwtdp.dat'
prop_azi = 90.0 # Propagation azimuth
isotherm_temp = 273.15+15. # Temperature of isothermal atmos in K
#----------------------------------------------------------------

df = pd.read_csv(filename,skiprows=None,names=['z','u','v','w','t','d','p'],sep=r"\s*")

azirad = prop_azi*np.pi/180.
alongpathwind = (df['u'].values*np.sin(azirad)) + (df['v'].values*np.cos(azirad))

abss = np.sqrt(402.8*df['t'].values) # Adiabatic Sound Speed
ess = alongpathwind+abss # Effective Sound Speed
isoss = np.sqrt(402.8*isotherm_temp)*np.ones(len(abss))

#-------------------------
# Print out to ASCII files
#-------------------------

atmos2deff_file = 'atmos_2d_eff_sound_sp.dat'
atmos2dadiab_file = 'atmos_2d_adiabatic_sound_sp.dat'
atmos2disotherm_file = 'atmos_2d_isothermal_sound_sp.dat'


fe2d = open(dirstr+atmos2deff_file, 'w')
for x in range(len(ess)):
    # Write out z in metres (*1000.)
    fe2d.write('{:6.0f} {:7.3f}\n'.format(df.iloc[x]['z']*1000.,ess[x]))
fe2d.close()

fa2d = open(dirstr+atmos2dadiab_file, 'w')
for x in range(len(abss)):
    # Write out z in metres (*1000.)
    fa2d.write('{:6.0f} {:7.3f}\n'.format(df.iloc[x]['z']*1000.,abss[x]))
fa2d.close()

fi2d = open(dirstr+atmos2disotherm_file, 'w')
for x in range(len(isoss)):
    # Write out z in metres (*1000.)
    fi2d.write('{:6.0f} {:7.3f}\n'.format(df.iloc[x]['z']*1000.,isoss[x]))
fi2d.close()

# Figure to show the three atmospheric profiles

para = {'axes.labelsize': 16, 'text.fontsize': 16, 'legend.fontsize': 14, 'xtick.labelsize': 14,'ytick.labelsize': 14, 'figure.subplot.left': 0.12, 'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11, 'figure.subplot.top': 0.97}
plt.rcParams.update(para)

fig = plt.figure(figsize=(6,8))
ax1 = fig.add_axes([0.1,0.1,0.75,0.75])

altv= df['z']*1000.

ax1.plot(ess,altv,'k-',label='Eff. Sd. Sp.')
ax1.plot(abss,altv,'k--',label='Adiab. Sd. Sp.')
ax1.plot(isoss,altv,'k:',label='Iso. Sd. Sp.')
ax1.set_xlabel('Sound Speed (m/s)')
ax1.set_ylabel('Altitude (m)')
ax1.legend(loc=4)
plt.savefig('example_atmospheres.png',bbox_inches='tight')
