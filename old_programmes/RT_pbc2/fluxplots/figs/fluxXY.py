import numpy as np
import matplotlib.pyplot as plt
from sys import exit
from scipy.integrate import simps

dirname = "../sim1/data/"
L = 10.
alpha = 0
d = 1.
v = 0.

x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
Nx = x.shape[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
Ny = y.shape[0]
X,Y = np.meshgrid(x,y)

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]




jx0 = np.loadtxt(dirname+"jx0.dat", delimiter=';')
jx1 = np.loadtxt(dirname+"jx1.dat", delimiter=';')
jy0 = np.loadtxt(dirname+"jy0.dat", delimiter=';')
jy1 = np.loadtxt(dirname+"jy1.dat", delimiter=';')


Jx0 = np.zeros((Nx,Ny))
Jy0 = np.zeros((Nx,Ny))
Jx1 = np.zeros((Nx,Ny))
Jy1 = np.zeros((Nx,Ny))



for xi in range(0,Nx):
	for yi in range(0,Ny):

		Jx0[xi][yi] = 0.5*( jx0[xi][yi] + jx0[xi+1][yi] )
		Jy0[xi][yi] = 0.5*( jy0[xi][yi] + jy0[xi][yi+1] )
		Jx1[xi][yi] = 0.5*( jx1[xi][yi] + jx1[xi+1][yi] )
		Jy1[xi][yi] = 0.5*( jy1[xi][yi] + jy1[xi][yi+1] )

J0 = np.sqrt(Jx0**2 + Jy0**2)
J1 = np.sqrt(Jx1**2 + Jy1**2)
plt.imshow(J1.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal')

nx = 3
ny = 1
X,Y = np.meshgrid(x[::nx],y[::ny])
Nx = X.shape[1]
Ny = Y.shape[0]

Jx0 = np.zeros((Nx,Ny))
Jy0 = np.zeros((Nx,Ny))
Jx1 = np.zeros((Nx,Ny))
Jy1 = np.zeros((Nx,Ny))


print(Nx,Ny)
th = 0.5
d = 1.
for xi in range(Nx):
    for yi in range(Ny):

        Jx0[xi][yi] = 0.5*( jx0[xi*nx][yi*ny] + jx0[xi*nx+1][yi*ny] )
        Jy0[xi][yi] = 0.5*( jy0[xi*nx][yi*ny] + jy0[xi*nx][yi*ny+1] )

        norm = np.sqrt(Jx0[xi][yi]**2 + Jy0[xi][yi]**2)
        
        if norm > th:
            Jx0[xi][yi] /= norm/d
            Jy0[xi][yi] /= norm/d
        else:
            Jx0[xi][yi] = 0
            Jy0[xi][yi] = 0

        Jx1[xi][yi] = 0.5*( jx1[xi*nx][yi*ny] + jx1[xi*nx+1][yi*ny] )
        Jy1[xi][yi] = 0.5*( jy1[xi*nx][yi*ny] + jy1[xi*nx][yi*ny+1] )

        norm = np.sqrt(Jx1[xi][yi]**2 + Jy1[xi][yi]**2)
        
        if norm > th:
            Jx1[xi][yi] /= norm/d
            Jy1[xi][yi] /= norm/d
        else:
            Jx1[xi][yi] = 0
            Jy1[xi][yi] = 0



scale = 30
hw = 10.
hl = 10.
plt.quiver(X,Y,Jx1.T,Jy1.T, color="red", scale=scale, headwidth=hw, headlength=hl)
plt.colorbar()
plt.gca().set_aspect("auto")

plt.ylim([-.75,.75])
plt.tight_layout()
plt.show()


