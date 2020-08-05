import numpy as np
import matplotlib.pyplot as plt
from sys import exit

dirname = "data/"
n = 0

x = np.loadtxt(dirname+"x.dat", delimiter=';')

r = np.loadtxt(dirname+"r{}.dat".format(n), delimiter=";")



plt.imshow(r)
plt.show()



