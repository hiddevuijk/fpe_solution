import numpy as np
import matplotlib.pyplot as plt
from sys import exit

#plt.imshow(data)
#plt.colorbar()

data = np.loadtxt("out.dat")
x = data[:,0]
rho = data[:,1]
print(sum(rho))

plt.plot(x,rho)
plt.show()



