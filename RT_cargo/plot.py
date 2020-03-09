import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

L = 10.
alpha = 0
d = 1.
v = 0.

r = np.loadtxt("r.dat", delimiter=';')
#s = np.loadtxt("s.dat", delimiter=';')
#print(sum(sum(r)))

#plt.subplot(1,2,1)
plt.imshow(r.T, origin='lower')
plt.colorbar()

#plt.subplot(1,2,2)
#plt.imshow(s.T, origin='lower')
#plt.colorbar()

plt.tight_layout()
plt.show()

