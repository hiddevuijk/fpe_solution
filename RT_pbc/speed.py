import numpy as np
import matplotlib.pyplot as plt
from ratio import mid, left
T = 1000

dirname = "data/"

L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]





m = []

for t in range(T):
	try:
		r = np.loadtxt(dirname+"r{}.dat".format(int(t) ), delimiter=';')
		ry = np.sum(r, axis = 0 )*dx
		rx = np.sum(r, axis = 1 )*dy
		m.append(mid(rx,dx))
	except:
		continue

t = np.loadtxt(dirname+"t.dat")
plt.plot(t,m)
plt.show()

