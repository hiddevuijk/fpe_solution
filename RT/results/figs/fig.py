import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from norm import norm_array

from chemtax import cn 
from functions import rho as rhoR_analytical
from functions import rho2 as rho2R_analytical

from get_data import get_data


q =  4.
kg = 1.
v0 = 10.
vp = 10.
x0 = 0
alpha = 40
D = 1.

name = "../q02_k1_v5_D1"

rho, R = get_data(name)

dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


def v(x, v0, vp, x0):
    return v0 + vp*np.sin(x*np.pi*2/L)

rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)
vlist = abs(vlist)
vavg = dR*np.sum(vlist)/(L)


rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L - 1, label="Numerical")

plt.plot(R,rhoA*L - 1, label="Analytical")
plt.plot(R,rho2A*L - 1, label="Analytical2")

plt.ylabel(r"$\frac{\rho(R)}{\rho_b}$", rotation=0, fontsize=10,labelpad=10)

plt.xlabel(r"$R$")
plt.legend()
plt.tight_layout()
plt.show()


