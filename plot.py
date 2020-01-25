import numpy as np
import matplotlib.pyplot as plt
from sys import exit

data = np.loadtxt("out.dat")

plt.imshow(data)
plt.colorbar()
plt.show()



