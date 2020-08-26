import numpy as np



def get_data(name):
    rho = np.loadtxt(name+"/r.dat", delimiter=';')
    R = np.loadtxt(name+"/x.dat")

    return rho, R

