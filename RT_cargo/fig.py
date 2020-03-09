import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt("x.dat")
dx = x[1] - x[0]
y = np.loadtxt("y.dat")
dy = y[1] - y[0]

r = np.loadtxt("r.dat", delimiter=';')
ry = np.sum(r, axis = 0 ) 
rx = np.sum(r, axis = 1 ) 
#s = np.mean(np.loadtxt("s.dat", delimiter=';'), axis = 1) 

norm = simps(ry,y)
ry /= norm
ryavg = simps(y*ry, y)

norm = simps(rx,x)
rx /= norm
rxavg = simps(x*rx, x)

plt.subplot(2,1,1)
plt.plot(y,ry, label="r")
plt.axvline(ryavg)
plt.legend()
plt.title(ryavg)

plt.subplot(2,1,2)
plt.plot(x,rx, label="r")
plt.axvline(rxavg)
plt.legend()
plt.title(rxavg)


plt.tight_layout()
plt.show()

