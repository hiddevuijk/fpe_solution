import numpy as np
import copy


def rho(v,dx, alpha, D, q):
    r = np.copy(v)**2
    r /= (1+q)*alpha*D
    r += 1
    r = r**(-(1-q)/2)
    N = np.sum(r)*dx 
    r /= N
    
    return r


    
