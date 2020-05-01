import numpy as np
import matplotlib.pyplot as plt
from sys import exit

L = 5.

v0 = 5.
vp = 2.

n = 4

N = 1000
X = np.linspace(-n*L/2,n*L/2,N)

def v(x):
	xx = x
	if(xx>L/2): xx = xx - L
	elif(xx<-L/2): xx = xx +  L
	if(xx>L/2): xx = xx - L
	elif(xx<-L/2): xx = xx +  L

	return v0-vp*abs(xx)

Y = np.asarray( [ v(xi) for xi in X ])


plt.plot(X,Y)
plt.axvline(L/2)
plt.axvline(-L/2)
plt.show()
