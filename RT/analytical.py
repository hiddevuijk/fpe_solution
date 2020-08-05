import numpy as np
from scipy.integrate import simps
from sys import exit


def rho(x,alpha, v, d, L):
    xi = 1./np.sqrt(v*v/(d*d) + alpha/d)

    a = alpha*d/v*v
    a *= np.exp( L/(2*xi) ) + np.exp( -L/(2*xi) )

    y = a + np.exp(-x/xi) +  np.exp(x/xi)

    norm = simps(y,x)
    y /= norm

    return y
    
