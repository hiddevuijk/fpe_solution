import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

dirname = "data/"
L = 10.
alpha = 0
d = 1.
v = 0.

r = np.loadtxt(dirname+"r.dat", delimiter=';')
x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
y = np.loadtxt(dirname+"y.dat")
dx = x[1] - x[0]
s = np.loadtxt(dirname+"s.dat", delimiter=';')

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]

plt.subplot(1,2,1)
#plt.imshow(r.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='equal')
plt.imshow(r.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
#plt.ylim([-2,2])
plt.title("density")
plt.xlabel("R")
plt.ylabel("r")
plt.colorbar(fraction=0.046)


plt.subplot(1,2,2)
#plt.imshow(s.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='equal')
plt.imshow(s.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
#plt.ylim([-2,2])
plt.title("otientation")
plt.xlabel("R")
plt.ylabel("r")
plt.colorbar(fraction=0.046)

plt.tight_layout()
plt.show()

