import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt("x.dat")
dx = x[1] - x[0]
y = np.loadtxt("y.dat")
dy = y[1] - y[0]

c = ['red', 'black', 'green', 'blue', 'orange', 'magenta', 'yellow', 'lightblue', 'pink']


for ni in range(1,10):
    r = np.loadtxt("r"+str(ni)+".dat", delimiter=';')
    ry = np.sum(r, axis = 0 ) 

    norm = simps(ry,y)
    ry /= norm
    ryavg = simps(y*ry, y)


    plt.plot(y,ry, label="r " + str(ni), color=c[ni%len(c)])
    plt.axvline(ryavg, color=c[ni%len(c)])


plt.legend()
plt.title(ryavg)

plt.tight_layout()
plt.show()

