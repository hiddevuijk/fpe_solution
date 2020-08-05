import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

dirname = "data2/"
L = 10.
alpha = 0
d = 1.
v = 0.

r = np.loadtxt(dirname+"r.dat", delimiter=';')
x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
Nx = x.shape[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
Ny = y.shape[0]
s = np.loadtxt(dirname+"s.dat", delimiter=';')

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]


rXY = np.zeros((Nx,Ny))
for xi in range(Nx):
	for yi in range(Ny):

		X = x[xi] + y[yi] + L
		Y = x[xi] - y[yi] + L
		
		Xi = int(X/(2*dx))
		Yi = int(Y/(2*dy))

		#if Xi >= 0 and Xi < Nx and Yi >= 0 and Yi < Ny:
		#	rXY[Xi][Yi] = r[xi][yi]
		rXY[Xi][Yi] = r[xi][yi]
		

plt.imshow(rXY.T,origin='lowerleft', interpolation='none', aspect='auto')
plt.colorbar()


plt.tight_layout()
plt.savefig("fig3.pdf")
plt.show()

