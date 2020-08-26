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

def rho2(v,dx, alpha, D, q, kg):
    r = np.copy(v)**2
    r /= (1+q)*alpha*D
    r += 1
    a = kg/alpha
    e = 1 - q*q/( q + (1+q)*a )
    r = r**(-e/2)
    N = np.sum(r)*dx 
    r /= N
    
    return r



