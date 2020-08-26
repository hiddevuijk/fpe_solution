import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps
from norm import norm_array

from chemtax import cn 
from functions import rho as rhoR_analytical
from functions import rho2 as rho2R_analytical

from get_data import get_data



name = "../q1.0488_k1_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]

rhoR = np.sum(rho,axis=1)
N = np.sum(rhoR)*dR
rhoR /=N

plt.plot(R,rhoR*L, label="k=1")



name = "../q1.1844_k4_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]

rhoR = np.sum(rho,axis=1)
N = np.sum(rhoR)*dR
rhoR /=N

plt.plot(R,rhoR*L, label="k=4")


name = "../q1.3483_k8_v5_D1"

rho, R = get_data(name)
dR = R[1] - R[0]
L = 2*R[-1]+dR
NR = rho.shape[0]

rhoR = np.sum(rho,axis=1)
N = np.sum(rhoR)*dR
rhoR /=N

plt.plot(R,rhoR*L, label="k=8")



plt.ylabel(r"$\frac{\rho(R)}{\rho_b}$", rotation=0, fontsize=20,labelpad=10)
plt.xlabel(r"$R$")

plt.title(r"$ q_0(k=1) = 1.0488 ~~ q_0(k=4)= 1.1844 ~~ q_0(k=8)= 1.3488$")

plt.legend()
plt.tight_layout()
plt.show()


