import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit

from chemtax import cn
from functions import rho as rhoR_analytical
from plot import plt_rho1d, plt_rho2d

fs = (8.27, 11.69)
fig, axes = plt.subplots(nrows=3,ncols=3, figsize=fs)

def v(x, v0, vp, x0,L):
    return abs(v0 + vp*np.sin(x*np.pi*2/L))
    #return abs(v0 - vp * abs(x-x0))



##############
#  q = 0.2
################
# R density
ax = axes[0,0]

simname = "../sim1/"

q = .2
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
D = 0.1

plt_rho1d(ax, q,v,alpha,D,simname)

# Rr density
ax = axes[0,1]

plt_rho2d(fig,ax, .25, simname)


# flux


##############
#  q = 1.
################
# R density
ax = axes[1,0]

simname = "../sim2/"

q = 1.
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
D = 0.1

plt_rho1d(ax, q,v,alpha,D,simname)

# Rr density
ax = axes[1,1]

plt_rho2d(fig,ax, .25, simname)


# flux



#######################

fig.tight_layout()
plt.show()

