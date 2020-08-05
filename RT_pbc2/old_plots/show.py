import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sys import exit
from scipy.integrate import simps

name = "fig2.pdf"
dirname = "data2/"
gamma = 5.
Gamma = 0.5
n = 3
L = 10.

x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
y = np.loadtxt(dirname+"y.dat")
dx = x[1] - x[0]

xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]

X,Y = np.meshgrid(x,y)

r = np.loadtxt(dirname+"r.dat", delimiter=';')
ry = np.sum(r, axis = 0 )*dx
rx = np.sum(r, axis = 1 )*dy

s = np.loadtxt(dirname+"s.dat", delimiter=';')

jx0 = np.loadtxt(dirname+"jx0.dat", delimiter=';')
jx1 = np.loadtxt(dirname+"jx1.dat", delimiter=';')
jy0 = np.loadtxt(dirname+"jy0.dat", delimiter=';')
jy1 = np.loadtxt(dirname+"jy1.dat", delimiter=';')

j0 = np.sqrt(jx0[:-1,:]**2 + jy0[:,:-1]**2)
j1 = np.sqrt(jx1[:-1,:]**2 + jy1[:,:-1]**2)

####################
#	figure
####################
figsize = (8.27,11.69)
fig, axes = plt.subplots(figsize=figsize, nrows=3, ncols=2)
fig.suptitle(r"$" + "\gamma = {} ~~~ \Gamma = {}".format(gamma, Gamma) + "$")


####################
# Density
####################

ax = axes[0,0]
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size='5%', pad = 0.05)

im = ax.imshow(r.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='equal')
#plt.ylim([-2,2])
ax.set_title("density")
ax.set_xlabel("R")
ax.set_ylabel("r")

fig.colorbar(im, cax=cax)

####################
# Orientation
####################

ax = axes[0,1]
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size='5%', pad = 0.05)

im = ax.imshow(s.T, extent=[xmin,xmax,ymin,ymax], interpolation='none', aspect='equal')
#plt.ylim([-2,2])
ax.set_title("Orientation")
ax.set_xlabel("R")
ax.set_ylabel("r")

fig.colorbar(im, cax=cax)

####################
# density x
####################

ax = axes[1,0]
ax.plot(x, rx)
ax.set_ylim([0,1.1*max(rx)])
ax.set_xlim([-L/2, L/2])
ax.set_ylabel(r'$\rho(R)$')
ax.set_xlabel(r'$R$')
ax.set_title("Density ABP")


####################
# density y
####################

ax = axes[1,1]
ax.plot(y, ry)
ax.set_ylim([0,1.1*max(ry)])
ax.set_xlim([-L/2, L/2])
ax.set_ylabel(r'$\rho(r)$')
ax.set_xlabel(r'$r$')
ax.set_title("Density BP")

####################
# density flux
####################
ax = axes[2,0]
divider  = make_axes_locatable(ax)
cax = divider.append_axes("right", size='5%', pad = 0.05)

ax.quiver(X[::n,::n],Y[::n,::n],jx0[:-1:n,::n].T,jy0[::n,:-1:n].T, color="red", headwidth=7., headlength=10.)
im = ax.imshow(j0.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal')
ax.set_title("density flux")
ax.set_xlabel("R")
ax.set_ylabel("r")


fig.colorbar(im, cax=cax)



####################
# density flux
####################
ax = axes[2,1]
divider  = make_axes_locatable(ax)
cax = divider.append_axes("right", size='5%', pad = 0.05)

ax.quiver(X[::n,::n],Y[::n,::n],jx1[:-1:n,::n].T,jy1[::n,:-1:n].T, color="red", headwidth=7., headlength=10.)
im = ax.imshow(j1.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal')
ax.set_title("orientation flux")
ax.set_xlabel("R")
ax.set_ylabel("r")


fig.colorbar(im, cax=cax)




####################
plt.tight_layout()
plt.subplots_adjust( top=0.94, bottom=0.05, left=0.095, right=0.945,hspace=0.272, wspace=0.24)
plt.savefig(name)
plt.show()

