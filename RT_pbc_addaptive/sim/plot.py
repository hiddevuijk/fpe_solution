import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from norm import norm_array

from small_friction import rho as rho_limit

q = 0.05

rho = np.loadtxt("data/r.dat", delimiter=';')
R = np.loadtxt("data/x.dat")
r = np.loadtxt("data/y.dat")
dR = R[1] - R[0]
dr = r[1] - r[0]
LR = 2*R[-1]+dR
Lr = 2*r[-1]+dr
NR = rho.shape[0]
Nr = rho.shape[1]


#plt.subplot(2,1,1)
# rho[Ri,ri]
plt.imshow(rho.T, origin='lowerleft', extent=[-LR/2,LR/2, -Lr/2, Lr/2], aspect='auto')
plt.colorbar()
plt.show()
exit()
L = 20
N = 1000
lp = 10

gamma = 1.
T = 1
k = 5
v0 = 0
vp = 5
x0 = 0
alpha = 20


xlist = np.linspace(-L/2,L/2,N)
d = xlist[1] - xlist[0]

def xy2R(x,y,q):
    return x/(1+q) + q*y/(1+q)
def xy2r(x, y,q):
    return x - y

def Rindex(x,y,q, dR):
    X = xy2R(x,y,q)
    return int( np.floor((X+LR/2)/dR) )
def rindex(x,y,q, dr):
    X = xy2r(x,y,q)
    return int( np.floor((X+Lr/2)/dr) )

rhoxy = np.zeros( (N,N) )

for i,x in enumerate(xlist):
    for j,y in enumerate(xlist):
        Ri = Rindex(x,y,q,dR)
        ri = rindex(x,y,q,dr)
        if Ri < 0: Ri += NR
        if Ri >= NR: Ri -= NR
        #if ri < 0: ri += Nr
        #if ri >= rho.shape[1]: ri -= Nr
        if Ri < NR and Ri >=0 and ri < Nr and ri>=0:
            rhoxy[i,j] = rho[Ri,ri]
        
            #if Ri == rho.shape[0] - 1:
            #    rhoxy[i,j] = 1
            #if Ri == 0:
            #    rhoxy[i,j] = 0.5
            #if ri == 0:
            #    rhoxy[i,j] = -0.5
            #if ri == rho.shape[1]-1:
            #    rhoxy[i,j] = -1
                
L/=2
#plt.subplot(2,1,1)
plt.imshow(rhoxy, origin='lowerleft', extent=[-L/2,L/2,-L/2,L/2], interpolation='none',aspect='auto')
plt.colorbar()
plt.show()


rhox = np.sum(rhoxy, axis=0)*d
rhoy = np.sum(rhoxy, axis=1)*d



xx, yy = norm_array(xlist,rhox,-L/2,L/2)
norm = simps(yy,xx)
yy /= norm
plt.plot(xx,yy, label=r"$\rho(x)$", color='red')


xx, yy = norm_array(xlist,rhoy,-L/2,L/2)
norm = simps(yy,xx)
yy /= norm
plt.plot(xx, yy, label=r"$\rho(y)$", color='blue')

rhox, rhoy, xlist = rho_limit(v0, vp, x0, alpha, T, gamma, q,k, L,lp,N)
xx, yy = norm_array(xlist,rhox,-L/2,L/2)
norm = simps(yy,xx)
yy /= norm
plt.plot(xx,yy, label=r"$\rho(x)$", color='red', ls=':')

xx, yy = norm_array(xlist,rhoy,-L/2,L/2)
norm = simps(yy,xx)
yy /= norm
plt.plot(xx,yy, label=r"$\rho(y)$", color='blue', ls=':')


plt.xlim([-L/2,L/2])
plt.legend()
plt.show()




