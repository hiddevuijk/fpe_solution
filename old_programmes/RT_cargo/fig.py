import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from sys import exit

dirname = "data/"

L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]

r = np.loadtxt(dirname+"r.dat", delimiter=';')
ry = np.sum(r, axis = 0 )*dx
rx = np.sum(r, axis = 1 )*dy
#s = np.mean(np.loadtxt("s.dat", delimiter=';'), axis = 1) 

norm = simps(ry,y)
#ry /= norm
ryavg = simps(y*ry, y)
#print(ryavg)

norm = simps(rx,x)
#rx /= norm
rxavg = simps(x*rx, x)

plt.subplot(2,1,1)
plt.plot(y,ry, label="r")
plt.axvline(ryavg)

plt.ylim([0,1.1*max(ry)])
plt.xlim([min(y),max(y)])

plt.legend()
plt.title(ryavg)


plt.subplot(2,1,2)
plt.plot(x,rx, label="r")
plt.axvline(rxavg)

plt.ylim([0,1.1*max(rx)])
plt.xlim([min(x),max(x)])
plt.legend()
plt.title(rxavg)


plt.tight_layout()
plt.savefig("fig.pdf")
plt.show()

