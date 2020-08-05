import numpy as np



def cn(x,rho,v, n=1):
    dx = x[1] - x[0]
    L = 2*x[-1] + dx

    rho_b = dx*np.sum(rho)/L
    vn = abs(v)
    vn = vn**n 
    vn_avg = dx*np.sum(vn)/L

    return dx*np.sum( ( (vn/vn_avg) - 1)*( (rho/rho_b) - 1) )/L
