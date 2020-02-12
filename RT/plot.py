import numpy as np
import matplotlib.pyplot as plt
from sys import exit

data = np.loadtxt("out.dat")

plt.subplot(1,2,1)
plt.ylim([0,1.])
plt.plot(data[:,0], label="rho")
plt.subplot(1,2,2)
plt.plot(data[:,1], label="sigma")
plt.show()
