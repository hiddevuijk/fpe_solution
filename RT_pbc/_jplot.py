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

jx0 = np.loadtxt(dirname+"jx0.dat", delimiter=';')
jx1 = np.loadtxt(dirname+"jx1.dat", delimiter=';')
jy0 = np.loadtxt(dirname+"jy0.dat", delimiter=';')
jy1 = np.loadtxt(dirname+"jy1.dat", delimiter=';')

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]

plt.subplot(2,2,1)
plt.imshow(jx0.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("jx0")
plt.colorbar()

plt.subplot(2,2,3)
plt.imshow(jx1.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("jx1")
plt.colorbar()

plt.subplot(2,2,2)
plt.imshow(jy0.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("jy0")
plt.colorbar()

plt.subplot(2,2,4)
plt.imshow(jy1.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.title("jy1")
plt.colorbar()

plt.tight_layout()
plt.show()

