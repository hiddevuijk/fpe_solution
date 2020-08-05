import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

dirname = "data/"
L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
y = np.loadtxt(dirname+"y.dat")
dx = x[1] - x[0]
X,Y = np.meshgrid(x,y)

jx0 = np.loadtxt(dirname+"jx0.dat", delimiter=';')
jx1 = np.loadtxt(dirname+"jx1.dat", delimiter=';')
jy0 = np.loadtxt(dirname+"jy0.dat", delimiter=';')
jy1 = np.loadtxt(dirname+"jy1.dat", delimiter=';')

j0 = np.sqrt(jx0[:-1,:]**2 + jy0[:,:-1]**2)
j1 = np.sqrt(jx1[:-1,:]**2 + jy1[:,:-1]**2)
'''
jx0 = np.loadtxt(dirname+"jx0.dat", delimiter=';').T
jx1 = np.loadtxt(dirname+"jx1.dat", delimiter=';').T
jy0 = np.loadtxt(dirname+"jy0.dat", delimiter=';').T
jy1 = np.loadtxt(dirname+"jy1.dat", delimiter=';').T

j0 = np.sqrt(jx0[:,:-1]**2 + jy0[:-1,:]**2)
j1 = np.sqrt(jx1[:,:-1]**2 + jy1[:-1,:]**2)
'''


xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]

#plt.subplot(1,2,1)
plt.quiver(X,Y,jx0[:-1,:].T,jy0[:,:-1].T, color="red", minlength=1.)
plt.imshow(j0.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal', origin='lower')
plt.colorbar()
'''
plt.subplot(1,2,2)
plt.quiver(X,Y,jx1[:-1,:],jy1[:,:-1], color="red", minlength=1.)
plt.imshow(j1, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='equal')
plt.colorbar()
'''

plt.tight_layout()
plt.show()


exit()
plt.subplot(1,2,1)
plt.imshow(j0.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("j0")
plt.colorbar()

plt.subplot(1,2,2)
plt.imshow(j1.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("j1")
plt.colorbar()

plt.tight_layout()
plt.show()

