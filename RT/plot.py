import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from analytical import rho

L = 5.
alpha = .1
d = .1
v = 1.


data = np.loadtxt("out.dat")

x = np.linspace(-L/2., L/2, 1000) 
y = rho(x,alpha, v, d,L)

norm = simps(data[:,1], data[:,0])

#plt.subplot(1,2,1)
plt.ylim([0,1.])

plt.plot(data[:,0],data[:,1]/norm, label="rho sim")
plt.plot(x,y)
plt.legend()

plt.show()
exit()
plt.subplot(1,2,2)
plt.plot(data[:,0], data[:,2], label="sigma")
plt.show()
