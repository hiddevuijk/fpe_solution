import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from norm import norm_array

from chemtax import cn 
from functions import rho as rhoR_analytical
from functions import rho2 as rho2R_analytical

from get_data import get_data


v0 = 5.
vp = 5.
x0 = 0
alpha = 40
D = 1.

plt.figure(figsize=(8.27,11.69) )

def v(x, v0, vp, x0):
    return v0 + vp*np.sin(x*np.pi*2/L)

##############################
# k = 1
##############################
kg=1
plt.subplot(3,1,1)
plt.title("k=1")

q=0.2
name = "../q02_k1_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='black',  label="q=0.2")
#plt.plot(R,rhoA*L , ls = ':', color='black' )
plt.plot(R,rho2A*L , ls = '--', color='black')

q=1
name = "../q1_k1_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='blue',  label="q=1.")
#plt.plot(R,rhoA*L , ls = ':', color='blue' )
plt.plot(R,rho2A*L , ls = '--', color='blue')


q=5
name = "../q5_k1_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='red',  label="q=5.")
#plt.plot(R,rhoA*L , ls = ':', color='black' )
plt.plot(R,rho2A*L , ls = '--', color='red')


plt.ylabel(r"$\frac{\rho(R)}{\rho_b}$", rotation=0, fontsize=10,labelpad=10)
plt.xlabel(r"$R$")
plt.xlim([-L/2,L/2])
plt.legend()

##############################
# k = 4
##############################
kg=4
plt.subplot(3,1,2)
plt.title("k=4")

q=0.2
name = "../q02_k4_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='black',  label="q=0.2")
#plt.plot(R,rhoA*L , ls = ':', color='black' )
plt.plot(R,rho2A*L , ls = '--', color='black')

q=1
name = "../q1_k4_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='blue',  label="q=1.")
#plt.plot(R,rhoA*L , ls = ':', color='blue' )
plt.plot(R,rho2A*L , ls = '--', color='blue')


q=5
name = "../q5_k4_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='red',  label="q=5.")
#plt.plot(R,rhoA*L , ls = ':', color='red' )
plt.plot(R,rho2A*L , ls = '--', color='red')

plt.ylabel(r"$\frac{\rho(R)}{\rho_b}$", rotation=0, fontsize=10,labelpad=10)
plt.xlabel(r"$R$")
plt.xlim([-L/2,L/2])

##############################
# k = 8
##############################
kg=8
plt.subplot(3,1,3)
plt.title("k=8")

q=0.2
name = "../q02_k8_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='black',  label="q=0.2")
#plt.plot(R,rhoA*L , ls = ':', color='black' )
plt.plot(R,rho2A*L , ls = '--', color='black')

q=1
name = "../q1_k8_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='blue',  label="q=1.")
#plt.plot(R,rhoA*L , ls = ':', color='blue' )
plt.plot(R,rho2A*L , ls = '--', color='blue')


q=5
name = "../q5_k8_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]


rhoR = np.sum(rho,axis=1)
rhoavg = dR*np.sum(rhoR)/L
N = np.sum(rhoR)*dR
rhoR /=N

vlist = v(R,v0,vp,x0)

rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)


plt.plot(R,rhoR*L ,color='red',  label="q=5.")
#plt.plot(R,rhoA*L , ls = ':', color='black' )
plt.plot(R,rho2A*L , ls = '--', color='red')

plt.ylabel(r"$\frac{\rho(R)}{\rho_b}$", rotation=0, fontsize=10,labelpad=10)
plt.xlabel(r"$R$")
plt.xlim([-L/2,L/2])

plt.tight_layout()
plt.show()


