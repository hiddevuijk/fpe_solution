import numpy as np
import matplotlib.pyplot as plt


dt = np.loadtxt("data/dt.dat", delimiter=';')
plt.plot(dt)
plt.show()
