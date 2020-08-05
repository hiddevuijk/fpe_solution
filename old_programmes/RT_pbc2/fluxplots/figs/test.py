import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit

from chemtax import cn
from functions import rho as rhoR_analytical
from plot import plt_rho1d, plt_rho2d, plt_sig1d, plt_sig2d, plt_j0, plt_j1


k = 8.
q = 0.2
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
D = 0.1
T = D

def v(x, v0, vp, x0,L):
    
    return v0 + vp*np.sin((x-x0)*np.pi*2/L)
    #return abs(v0 + vp*np.sin(x*np.pi*2/L))
    #return abs(v0 - vp * abs(x-x0))

simname = "../sim1/"
rho = np.loadtxt(simname+"data/r.dat", delimiter=';')
R = np.loadtxt(simname+"data/x.dat")
r = np.loadtxt(simname+"data/y.dat")
dR = R[1] - R[0]
dr = r[1] - r[0]
L = R[-1] - R[0] + dR
print(L)
rhoR = np.sum(rho,axis=1)*dr

vlist = v(R,v0,vp,x0,L)
c1 = cn(R,rhoR,vlist,1)
c2 = cn(R,rhoR,vlist,2)

print(c1,c2)
plt.plot(R,vlist)
plt.show()


