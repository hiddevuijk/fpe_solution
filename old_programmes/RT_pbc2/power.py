import numpy as np
import matplotlib.pyplot as plt
from sys import exit

q = 5.
k = 5.
gamma = 1
T = 1

rho = np.loadtxt("data/r.dat", delimiter=';')
R = np.loadtxt("data/x.dat")
r = np.loadtxt("data/y.dat")
dR = R[1] - R[0]
dr = r[1] - r[0]

N = sum(sum(rho))*dR*dr
rho /= N


w = -k*T/(q*gamma)
for Ri in range(rho.shape[0]):
    for ri in range(rho.shape[1]):
        dw = rho[Ri,ri] * r[ri]**2
        dw *= dr*dR
        dw *= k*k/(q*gamma)
        w += dw

print(w)


