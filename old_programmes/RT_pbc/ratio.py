import numpy as np
from sys import exit

def mid(r,d):
	norm = np.sum(r)*d

	N = r.shape[0]
	if N%2 ==0 :
		m = np.sum( r[ int(N/4) : int(3*N/4) ] )*d
		
	return m/norm

def left(r,d):
	norm = np.sum(r)*d
	N = r.shape[0]
	if N%2 == 0:
		m = np.sum( r[0: int(N/4)] )*d
		m += np.sum( r[ int(3*N/4):] )*d
	return m/norm



