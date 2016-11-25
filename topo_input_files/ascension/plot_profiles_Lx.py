#! /usr/local/python27/bin/python


# Short script to read in the GMT output data and plot 
# four profiles to the different stations.


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_profile(filename):
	data = pd.io.parsers.read_csv(filename, sep=r'\s*',names=['lon','lat','dist','alt'])
	return data
	
para = {'axes.labelsize': 18, 'text.fontsize': 18, 'legend.fontsize': 13, 'xtick.labelsize': 16,'ytick.labelsize': 16, 'figure.subplot.left': 0.12, 'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11, 'figure.subplot.top': 0.97}
plt.rcParams.update(para)

L1data = read_profile('L1.xydz')
L2data = read_profile('L2.xydz')
L3data = read_profile('L3.xydz')
L4data = read_profile('L4.xydz')

fig = plt.figure(figsize=(10,10))
ax3 = fig.add_axes([0.1,0.5,0.8,0.15])
ax1 = fig.add_axes([0.1,0.1,0.8,0.15],sharex=ax3,sharey=ax3)
ax2 = fig.add_axes([0.1,0.3,0.8,0.15],sharex=ax3,sharey=ax3)
ax4 = fig.add_axes([0.1,0.7,0.8,0.15],sharex=ax3,sharey=ax3)

ax1.plot(L1data['dist'],L1data['alt'],'k-')
ax1.plot(L1data['dist'].iloc[0],L1data['alt'].iloc[0],'r*',ms=12,clip_on=False)
ax1.plot(L1data['dist'].iloc[-1],L1data['alt'].iloc[-1],'b^',ms=12,clip_on=False)
ax1.text(L1data['dist'].iloc[-1]+0.1,L1data['alt'].iloc[-1],'L1',fontsize=20)

ax2.plot(L2data['dist'],L2data['alt'],'k-')
ax2.plot(L2data['dist'].iloc[0],L2data['alt'].iloc[0],'r*',ms=12,clip_on=False)
ax2.plot(L2data['dist'].iloc[-1],L2data['alt'].iloc[-1],'b^',ms=12,clip_on=False)
ax2.text(L2data['dist'].iloc[-1]+0.1,L2data['alt'].iloc[-1],'L2',fontsize=20)

ax3.plot(L3data['dist'],L3data['alt'],'k-')
ax3.plot(L3data['dist'].iloc[0],L3data['alt'].iloc[0],'r*',ms=12,clip_on=False)
ax3.plot(L3data['dist'].iloc[-1],L3data['alt'].iloc[-1],'b^',ms=12,clip_on=False)
ax3.text(L3data['dist'].iloc[-1]+0.1,L3data['alt'].iloc[-1],'L3',fontsize=20)

ax4.plot(L4data['dist'],L4data['alt'],'k-')
ax4.plot(L4data['dist'].iloc[0],L4data['alt'].iloc[0],'r*',ms=12,clip_on=False)
ax4.plot(L4data['dist'].iloc[-1],L4data['alt'].iloc[-1],'b^',ms=12,clip_on=False)
ax4.text(L4data['dist'].iloc[-1]+0.1,L4data['alt'].iloc[-1],'L4',fontsize=20)

ax1.set_xlabel('Along Path Distance (km)')
ax2.set_ylabel('Altitude (m a.s.l.)')


#-------------------------------------------
plt.savefig('ascension_profiles.png',bbox_inches='tight')


#------------------------------------------
# Now save to ASCII files, ensuring topography is in metres (not kilometres)

dirstr = '/home/dgreen/local/ascension/topographic_set/'
rzL1file = 'asc_L1_2d.dat'
rzL2file = 'asc_L2_2d.dat'
rzL3file = 'asc_L3_2d.dat'
rzL4file = 'asc_L4_2d.dat'

fL1 = open(dirstr+rzL1file,'w')
for x in range(len(L1data['dist'])):
	fL1.write('{:8.1f} {:7.3f}\n'.format(L1data.iloc[x]['dist']*1000.,L1data.iloc[x]['alt']))
fL1.close()

fL2 = open(dirstr+rzL2file,'w')
for x in range(len(L2data['dist'])):
	fL2.write('{:8.1f} {:7.3f}\n'.format(L2data.iloc[x]['dist']*1000.,L2data.iloc[x]['alt']))
fL2.close()

fL3 = open(dirstr+rzL3file,'w')
for x in range(len(L3data['dist'])):
	fL3.write('{:8.1f} {:7.3f}\n'.format(L3data.iloc[x]['dist']*1000.,L3data.iloc[x]['alt']))
fL3.close()

fL4 = open(dirstr+rzL4file,'w')
for x in range(len(L4data['dist'])):
	fL4.write('{:8.1f} {:7.3f}\n'.format(L4data.iloc[x]['dist']*1000.,L4data.iloc[x]['alt']))
fL4.close()






