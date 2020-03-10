import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

L = 10.
alpha = 0
d = 1.
v = 0.

r = np.loadtxt("r.dat", delimiter=';')
x = np.loadtxt("x.dat")
y = np.loadtxt("y.dat")
#s = np.loadtxt("s.dat", delimiter=';')
#print(sum(sum(r)))

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]

#plt.subplot(1,2,1)
plt.imshow(r.T, origin='lower', extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='auto')
plt.colorbar()

#plt.subplot(1,2,2)
#plt.imshow(s.T, origin='lower')
#plt.colorbar()

plt.tight_layout()
plt.show()

