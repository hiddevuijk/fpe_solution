import numpy as np
import matplotlib.pyplot as plt
from sys import exit

dirname = "data/"

aList = [1/10.,  1/6., 2/6., 3/6., 4/6., 5/6.,1, 6/5.,6/4.,6/3., 6/2.,6., 10.]
#aList = [1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1,2.2,2.3,2.4,2.5]

x = np.loadtxt(dirname+"x.dat", delimiter=';')
y = np.loadtxt(dirname+"y.dat", delimiter=';')
dy = y[1] - y[0]

for ai, a in enumerate(aList):
	try:
		r = np.loadtxt(dirname+"rx{}.dat".format(ai), delimiter=";")
		plt.plot(x,r,label=a)
	except:
		continue
plt.legend()
plt.tight_layout()
plt.show()



