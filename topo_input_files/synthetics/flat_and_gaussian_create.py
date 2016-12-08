#! /usr/local/python27/bin/python

# For Codor -
# Setting up some 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy.integrate as integrate
import random as rdm

def gaussprof(A,hw,x,x0):
    
      # A function based on the matlab code of 
      # Codor Khodr (Uni. of Bristol) to generate the
      # a Gaussian hill topography profile
    
      # Input:
	# A = maximum altitude of the Gaussian function
	# hw = half-width (standard deviation) of Gaussian
	# x = input abssica
	# x0 = center of the Gaussian Hill
	
      # Output:
      # f -> the Gaussian hill function, i.e., f(x)
      
 
	# Generate the Gaussian
	
	n = x.shape[0]
	f = np.zeros(n)
	
      # Codor's function parameterised in terms of theta_max. Why?
	#A = np.tan(theta_max*np.pi/180.0)*d*np.exp(0.5)/np.sqrt(2)
	#print 'A = ',A
	
	for i in range(n):
		f[i] = A*np.exp(-1*((x[i]-x0)**2)/(2*(hw**2)))
		
	return f

def gaussprof3D(A,hwx,hwy,x,y,x0,y0):
    
    # Create a 3D Gaussian Hill
    # A function based on the matlab code of 
    # Codor Khodr (Uni. of Bristol) to generate the
    # a Gaussian hill topography profile
    
    # Input:
    # A = maximum height of the Gaussian Hill
    # hwx = half-width (standard deviation) of Gaussian in x-dirn
    # hwy = half-width (standard deviation) of Gaussian in y-dirn
    # x = input abssica
    # x0 = center of the Gaussian Hill
	
    # Output:
    # f -> the Gaussian hill function, i.e., f(x)

    # Create the mesh and calculate the Gaussian function
    xv,yv = np.meshgrid(x,y)
    # Flattening (inefficient programming step?)
    xvflat = xv.flatten()
    yvflat = yv.flatten()
    
    n = len(xvflat)
    f = np.zeros(n)
    
    for i in range(n):
        term1 = ((xvflat[i]-x0)**2)/(2*(hwx**2))
        term2 = ((yvflat[i]-y0)**2)/(2*(hwy**2))
        f[i]  = A*np.exp(-1*(term1+term2))
    
    return xvflat,yvflat,f

#----- Main Prog


#---------------------------------------------------------------
# Setting up the Gaussian, as per Codor's email of 10-March-2016
#
# Calculating path length using the standard formula
# p.29 of Riley, Hobson and Bence
#---------------------------------------------------------------
para = {'axes.labelsize': 16, 'text.fontsize': 16, 'legend.fontsize': 14, 'xtick.labelsize': 14,'ytick.labelsize': 14, 'figure.subplot.left': 0.12, 'figure.subplot.right': 0.98, 'figure.subplot.bottom': 0.11, 'figure.subplot.top': 0.97}
plt.rcParams.update(para)


#-----------
# 2D CASES
#-----------

L = 6000.
A = 200.
dx = 90.
x0 = L/2.0
hw = 3.0*A


x_regular = np.arange(0,L,dx)

# Calculating the regular topography function.

f_regular = gaussprof(A,hw,x_regular,x0)

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_axes([0.2,0.2,0.7,0.7])

ax1.plot(x_regular,f_regular,'ko')

plt.xlabel('x (m)')
plt.ylabel('y (m)')


# Calculating the irregularly sampled topography function

x_irregular = np.sort(rdm.sample(range(0,int(L)),len(x_regular)))

f_irregular = gaussprof(A,hw,x_irregular,x0)

ax1.plot(x_irregular,f_irregular,'r*')
plt.show()

#-----------------------------------------------------------------
# Now to output to ASCII files
#-----------------------------------------------------------------
dirpath = '/Users/dgreen/Documents/Work/4codor/infratopo/topo_input_files/synthetics/'

regfile = 'gauss_'+str(int(A))+'m_hill_short.dat'
irregfile = 'gauss_'+str(int(A))+'m_hill_irreg_short.dat'
flatfile = 'flat_topo_short.dat'

fr = open(dirpath+regfile, 'w')
for x in range(len(f_regular)):
    fr.write('{:4.0f} {:5.1f}\n'.format(x_regular[x],f_regular[x]))
fr.close()

fi = open(dirpath+irregfile, 'w')
for x in range(len(f_irregular)):
    fi.write('{:4.0f} {:5.1f}\n'.format(x_irregular[x],f_irregular[x]))
fi.close()

ff = open(dirpath+flatfile, 'w')
for x in range(len(f_regular)):
    ff.write('{:4.0f} {:5.1f}\n'.format(x_regular[x],0.0))
ff.close()


#-----------------------------------------------------------------
# Now for the long-range cases (500km)
#-----------------------------------------------------------------

L = 550000.
A = 3000.
dx = 1000.
x0 = 250000.
hw = 60000.


x_regular = np.arange(0,L,dx)

# Calculating the regular topography function.

f_regular = gaussprof(A,hw,x_regular,x0)

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_axes([0.2,0.2,0.7,0.7])

ax1.plot(x_regular,f_regular,'ko')

plt.xlabel('x (m)')
plt.ylabel('y (m)')

dirpath = '/Users/dgreen/Documents/Work/4codor/infratopo/topo_input_files/synthetics/'

longfile = 'gauss_'+str(int(A))+'m_hill_long.dat'
fr = open(dirpath+longfile, 'w')
for x in range(len(f_regular)):
    fr.write('{:6.0f} {:5.1f}\n'.format(x_regular[x],f_regular[x]))
fr.close()


longflat = 'flat_topo_long.dat'
fr = open(dirpath+longflat, 'w')
for x in range(len(f_regular)):
    fr.write('{:6.0f} {:5.1f}\n'.format(x_regular[x],0.0))
fr.close()

#-----------
# 3D CASES
#-----------

L = 6000.
A = 200.
dx = 90.
dy = dx
x0 = L/2.0
y0 = L/2.0
hwx = 3.0*A
hwy = 3.0*A

x_regular = np.arange(0,L,dx)
y_regular = np.arange(0,L,dy)

xv3d,yv3d,f3d = gaussprof3D(A,hwx,hwy,x_regular,y_regular,x0,y0)

# For plottting (not ascii output) need to convert above to 2D arrays

cols = np.unique(xv3d).shape[0]
X = xv3d.reshape(-1, cols)
Y = yv3d.reshape(-1, cols)
Z = f3d.reshape(-1, cols)


fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_axes([0.25, 0.25, 0.75, 0.75], projection='3d')
ax.set_zlabel('Alt (m)')
ax.set_ylabel('y (m)')
ax.set_xlabel('x (m)')
surf = ax.plot_surface(X,Y,Z,rstride=2,cstride=2,cmap=cm.coolwarm,linewidth=0, antialiased=False)
#fig.colorbar(surf, shrink=0.4, aspect=10)
ax.dist=11

plt.show()
fig.savefig('gauss_3D_example_crude.png',bbox_inches='tight')

# Output to ASCII
regfile3d = 'gauss3d_'+str(int(A))+'m_hill_short.dat'
flatfile3d = 'flat_topo_3d_short.dat'

fr3d = open(dirpath+regfile3d, 'w')
for x in range(len(f3d)):
    fr3d.write('{:4.0f} {:4.0f} {:5.1f}\n'.format(xv3d[x],yv3d[x],f3d[x]))
fr3d.close()

fl3d = open(dirpath+flatfile3d, 'w')
for x in range(len(f3d)):
    fl3d.write('{:4.0f} {:4.0f} {:5.1f}\n'.format(xv3d[x],yv3d[x],0.0))
fl3d.close()