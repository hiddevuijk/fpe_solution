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


Nx = x.shape[0]

n4 = int(Nx/4)


wl = simps( rx[:n4], x[:n4])
wr = simps( rx[n4-1:2*n4], x[n4-1:2*n4])
wr1 = simps( rx[2*n4-1:3*n4], x[2*n4-1:3*n4] )
wl1 = simps( rx[3*n4-1:], x[3*n4-1:] )

print(wl,wr)
print(wl1,wr1)
print(wl+wr+wr1+wl1)
print( simps( rx,x) )



