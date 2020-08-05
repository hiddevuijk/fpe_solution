import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from sys import exit

from chemtax import cn
from functions import rho as rhoR_analytical
from plot import plt_rho1d, plt_rho2d, plt_sig1d, plt_sig2d, plt_j0, plt_j1

fs = (8.27, 11.69)
fig, axes = plt.subplots(nrows=3,ncols=2, figsize=fs)

def v(x, v0, vp, x0,L):
    return abs(v0 + vp*np.sin(x*np.pi*2/L))
    #return abs(v0 - vp * abs(x-x0))



##############
################

simname = "../sim6k0/"
savename = "../../../../figures/q5_T1_k0"
q = 5.
D = 1.
k = 0.
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
T = D

# R density
ax = axes[0,0]
plt_rho1d(ax, q,v,alpha,D,simname)

# R sig
ax = axes[0,1]
plt_sig1d(ax,simname)

ry = .5
# Rr density
ax = axes[1,0]
plt_rho2d(fig,ax, ry, simname)

# Rr sig
ax = axes[1,1]
plt_sig2d(fig,ax, ry, simname)



# flux j0
ax = axes[2,0]
plt_j0(fig, ax, ry, 0.025, 4,5, simname)


# flux j1
ax = axes[2,1]
plt_j1(fig, ax, ry, 0.05, 4,5, simname)

########################
rho = np.loadtxt(simname+"data/r.dat", delimiter=';')
R = np.loadtxt(simname+"data/x.dat")
r = np.loadtxt(simname+"data/y.dat")
dR = R[1] - R[0]
dr = r[1] - r[0]
L = R[-1] - R[0] + dR
rhoR = np.sum(rho,axis=1)*dr



vlist = v(R,v0,vp,x0,L)
c1 = cn(R,rhoR,vlist,1)
c2 = cn(R,rhoR,vlist,2)
#######################
fig.suptitle(r"$q={:1.1f}~ T={:1.1f} k={:1.1f} ~~~~~ c_1={:1.3f}~c_2={:1.3f}$".format(q,D,k,c1,c2) )
fig.tight_layout()
fig.subplots_adjust(top=0.95)


plt.savefig(savename+".pdf")
plt.show()


