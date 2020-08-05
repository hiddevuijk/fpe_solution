import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import simps
from sys import exit

from chemtax import cn
from functions import rho as rhoR_analytical

simname = "../sim1/"

q = .2
v0 = 10.
vp = 10.
x0 = 0
alpha = 100
D = 0.1

def v(x, v0, vp, x0,L):
    return abs(v0 + vp*np.sin(x*np.pi*2/L))
    #return abs(v0 - vp * abs(x-x0))



def plt_rho1d(ax,q,v,alpha,D,simname):
    rho = np.loadtxt(simname+"data/r.dat", delimiter=';')
    R = np.loadtxt(simname+"data/x.dat")
    r = np.loadtxt(simname+"data/y.dat")
    dR = R[1] - R[0]
    dr = r[1] - r[0]
    L = 2*R[-1]+dR
    NR = rho.shape[0]
    Nr = rho.shape[1]

    rhoR = np.sum(rho,axis=1)
    rhoavg = dR*np.sum(rhoR)/L
    N = np.sum(rhoR)*dR
    rhoR /=N


    vlist = v(R,v0,vp,x0,L)

    c1 = cn(R,rhoR,vlist,1)
    c2 = cn(R,rhoR,vlist,2)


    rhoA = rhoR_analytical(vlist, dR,alpha,D,q)

    ax.plot(R,rhoR*L - 1, label="Numerical")
    ax.plot(R,rhoA*L - 1, label="Analytical")
    ax.set_ylabel(r"$\frac{\rho(R)}{\rho_b} - 1$", rotation=0, fontsize=10,labelpad=10)
    ax.set_xlabel(r"$R$")
    #ax.set_title( r"$q = {:1.1f} ~ D={:1.1f} ~~~~~ c_1 = {:1.3f}  ~ c_2={:1.3f} $".format(q,D, c1, c2))
    #ax.title( r"$q = {:1.2f} ~~ D={:1.1f}  $".format(q,D))
    ax.legend()

def plt_sig1d(ax,simname):
    sig = np.loadtxt(simname+"data/s.dat", delimiter=';')
    R = np.loadtxt(simname+"data/x.dat")
    r = np.loadtxt(simname+"data/y.dat")
    dR = R[1] - R[0]
    dr = r[1] - r[0]
    L = 2*R[-1]+dR

    sigR = np.sum(sig,axis=1)*dR

    ax.plot(R,sigR , label="Numerical")
    ax.set_ylabel(r"$\sigma(R)$", rotation=0, fontsize=10,labelpad=10)
    ax.set_xlabel(r"$R$")
 
def plt_rho2d(fig, ax, rlim , simname):
    rho = np.loadtxt(simname+"data/r.dat", delimiter=';')
    R = np.loadtxt(simname+"data/x.dat")
    dR = R[1] - R[0]
    LR = R[-1] - R[0] + dR
    r = np.loadtxt(simname+"data/y.dat")
    dr = r[1] - r[0]
    Lr = r[-1] - r[0] + dr
    im = ax.imshow(rho.T, origin='lowerleft', extent=[-LR/2,LR/2, -Lr/2, Lr/2])
    ax.set_ylim(-rlim,rlim)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$r$")
    ax.set_aspect('auto')
    ax.set_title(r"$\rho(R,r)$")
 
def plt_sig2d(fig, ax, rlim , simname):
    sig = np.loadtxt(simname+"data/s.dat", delimiter=';')
    R = np.loadtxt(simname+"data/x.dat")
    dR = R[1] - R[0]
    LR = R[-1] - R[0] + dR
    r = np.loadtxt(simname+"data/y.dat")
    dr = r[1] - r[0]
    Lr = r[-1] - r[0] + dr
    im = ax.imshow(sig.T, origin='lowerleft', extent=[-LR/2,LR/2, -Lr/2, Lr/2])
    ax.set_ylim(-rlim,rlim)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$r$")
    ax.set_aspect('auto')
    ax.set_title(r"$\sigma(R,r)$")

def plt_j0(fig, ax, rlim, th, nx,ny, simname):

    x = np.loadtxt(simname+"data/x.dat")
    dx = x[1] - x[0]
    Nx = x.shape[0]
    y = np.loadtxt(simname+"data/y.dat")
    dy = y[1] - y[0]
    Ny = y.shape[0]
    X,Y = np.meshgrid(x,y)

    xmin = x[0]
    xmax = x[-1]
    ymin = y[0]
    ymax = y[-1]




    jx0 = np.loadtxt(simname+"data/jx0.dat", delimiter=';')
    jy0 = np.loadtxt(simname+"data/jy0.dat", delimiter=';')


    Jx0 = np.zeros((Nx,Ny))
    Jy0 = np.zeros((Nx,Ny))

    for xi in range(0,Nx):
        for yi in range(0,Ny):

            Jx0[xi][yi] = 0.5*( jx0[xi][yi] + jx0[xi+1][yi] )
            Jy0[xi][yi] = 0.5*( jy0[xi][yi] + jy0[xi][yi+1] )

    J0 = np.sqrt(Jx0**2 + Jy0**2)
    im = ax.imshow(J0.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')


    X,Y = np.meshgrid(x[::nx],y[::ny])
    Nx = X.shape[1]
    Ny = Y.shape[0]

    Jx0 = np.zeros((Nx,Ny))
    Jy0 = np.zeros((Nx,Ny))

    for xi in range(Nx):
        for yi in range(Ny):

            Jx0[xi][yi] = 0.5*( jx0[xi*nx][yi*ny] + jx0[xi*nx+1][yi*ny] )
            Jy0[xi][yi] = 0.5*( jy0[xi*nx][yi*ny] + jy0[xi*nx][yi*ny+1] )

            norm = np.sqrt(Jx0[xi][yi]**2 + Jy0[xi][yi]**2)
            
            if norm > th:
                Jx0[xi][yi] /= norm
                Jy0[xi][yi] /= norm
            else:
                Jx0[xi][yi] = 0
                Jy0[xi][yi] = 0

    scale = 30
    hw = 10.
    hl = 10.
    ax.quiver(X,Y,Jx0.T,Jy0.T, color="red", scale=scale, headwidth=hw, headlength=hl)

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$r$")
    ax.set_ylim([-rlim,rlim])
    ax.set_aspect("auto")
    ax.set_title(r"$|J^{(0)}|$")
 
def plt_j1(fig, ax, rlim, th, nx,ny, simname):

    x = np.loadtxt(simname+"data/x.dat")
    dx = x[1] - x[0]
    Nx = x.shape[0]
    y = np.loadtxt(simname+"data/y.dat")
    dy = y[1] - y[0]
    Ny = y.shape[0]
    X,Y = np.meshgrid(x,y)

    xmin = x[0]
    xmax = x[-1]
    ymin = y[0]
    ymax = y[-1]

    jx1 = np.loadtxt(simname+"data/jx1.dat", delimiter=';')
    jy1 = np.loadtxt(simname+"data/jy1.dat", delimiter=';')

    Jx1 = np.zeros((Nx,Ny))
    Jy1 = np.zeros((Nx,Ny))

    for xi in range(0,Nx):
        for yi in range(0,Ny):

            Jx1[xi][yi] = 0.5*( jx1[xi][yi] + jx1[xi+1][yi] )
            Jy1[xi][yi] = 0.5*( jy1[xi][yi] + jy1[xi][yi+1] )

    J1 = np.sqrt(Jx1**2 + Jy1**2)
    im = ax.imshow(J1.T, extent=[xmin,xmax,ymin,ymax],interpolation='none', aspect='equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')


    X,Y = np.meshgrid(x[::nx],y[::ny])
    Nx = X.shape[1]
    Ny = Y.shape[0]

    Jx1 = np.zeros((Nx,Ny))
    Jy1 = np.zeros((Nx,Ny))

    for xi in range(Nx):
        for yi in range(Ny):

            Jx1[xi][yi] = 0.5*( jx1[xi*nx][yi*ny] + jx1[xi*nx+1][yi*ny] )
            Jy1[xi][yi] = 0.5*( jy1[xi*nx][yi*ny] + jy1[xi*nx][yi*ny+1] )

            norm = np.sqrt(Jx1[xi][yi]**2 + Jy1[xi][yi]**2)
            
            if norm > th:
                Jx1[xi][yi] /= norm
                Jy1[xi][yi] /= norm
            else:
                Jx1[xi][yi] = 0
                Jy1[xi][yi] = 0



    scale = 30
    hw = 10.
    hl = 10.
    ax.quiver(X,Y,Jx1.T,Jy1.T, color="red", scale=scale, headwidth=hw, headlength=hl)

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$r$")
    ax.set_aspect("auto")
    ax.set_ylim([-rlim,rlim])
    ax.set_title(r"$|J^{(1)}|$")

#fs = (8.27, 11.69)
#fig, axes = plt.subplots(nrows=3,ncols=3, figsize=fs)
#
#
#plt1(axes[0,0], q,v,alpha,D,simname)
#fig.tight_layout()
#
#plt.show()
#
