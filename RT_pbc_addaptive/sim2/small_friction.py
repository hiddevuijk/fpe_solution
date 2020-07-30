import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps, quad
from sys import exit

'''
q = 0.05
gamma = 1.
T = 1
k = 5
v0 = 5
vp = 5
x0 = 0
alpha = 20
L = 10
lp = 10
N = 1000
'''

def rho( v0,vp,x0, alpha, T, gamma, q, k, L, lp, N):
    def v(x):
        if x > lp/2: x-= lp
        if x < -lp/2: x+= lp
        #return v0*np.sin( np.pi*(x-x0)/L)
        return v0 - vp*abs(x - x0)

    def U(x,y):
        return 0.5*k*(x-y)**2


    def rho0(x):
        A = v(x)**2
        A /= alpha
        A *= gamma/T
        A *= 1-q
        return 1./np.sqrt(1+A)

    def rhoY(y):
        def I(x,y):
            return np.exp(-U(x,y)/T)*rho0(x)
        d = 5
        A = quad( lambda x: I(x,y), -L/2 - d,L/2 + d) 
        return A[0]


    x = np.linspace(-lp/2,lp/2,N)
    rx = [ rho0(xi) for xi in x ]
    norm = simps(rx,x)
    rx /= norm

    y = np.linspace(-lp/2,lp/2,N)
    ry =  [ rhoY(yi) for yi in y ]
    norm = simps(ry,y)
    ry /= norm
    return rx, ry, x


'''
plt.plot(x,rx, label="rho0")
plt.plot(y,ry, label="rhoY")

plt.legend()
plt.show()
'''
