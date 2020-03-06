import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

L = 5.
alpha = .1
d = 1.
v = 1.

r = np.loadtxt("r.dat", delimiter=';')

plt.imshow(r)
plt.colorbar()
plt.show()

