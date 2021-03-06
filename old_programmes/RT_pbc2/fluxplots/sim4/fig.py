import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

from chemtax import cn 
from functions import rho as rhoR_analytical


q = .2
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
D = 1.

rho = np.loadtxt("data/r.dat", delimiter=';')
R = np.loadtxt("data/x.dat")
r = np.loadtxt("data/y.dat")
dR = R[1] - R[0]
dr = r[1] - r[0]
L = 2*R[-1]+dR
NR = rho.shape[0]
Nr = rho.shape[1]

def v(x, v0, vp, x0):
    return v0 + vp*np.sin(x*np.pi*2/L)
    #return abs(v0 - vp * abs(x-x0))

rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)
vlist = abs(vlist)
vavg = dR*np.sum(vlist)/(L)


rhoA = rhoR_analytical(vlist, dR,alpha,D,q)

c = dR * np.sum( (vlist-vavg)*(rhoR -rhoavg) ) /( vavg*rhoavg* L)

c1 = cn(R,rhoR,abs(vlist),1)
c2 = cn(R,rhoR,vlist,2)
print(c1,c2)

plt.plot(R,rhoR*L - 1, label="Numerical")
plt.plot(R,rhoA*L - 1, label="Analytical")
plt.ylabel(r"$\frac{\rho(R)}{\rho_b} - 1$", rotation=0, fontsize=10,labelpad=10)
plt.xlabel(r"$R$")
#plt.title( r"$q = {:1.2f} ~~ D={:1.1f} ~~ c_1 = {:1.5f}  ~~ c_2={:1.5f} $".format(q,D, c1, c2))
plt.title( r"$q = {:1.2f} ~~ D={:1.1f}  $".format(q,D))
plt.legend()
plt.tight_layout()
plt.show()


