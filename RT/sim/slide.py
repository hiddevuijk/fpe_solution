import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from scipy.integrate import simps
from sys import exit

from functions import rho as rhoR_analytical
from functions import rho2 as rho2R_analytical

dirname = "data/"
Nfig = 105

fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.2)



x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
Lx = 2*x[-1]+dx
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]

#r = np.loadtxt(dirname+"r0.dat", delimiter=';')
#ry = np.sum(r, axis = 0 )*dx
#rx = np.sum(r, axis = 1 )*dy
rx = np.loadtxt(dirname+"rx0.dat", delimiter=';')
ry = np.loadtxt(dirname+"ry0.dat", delimiter=';')

norm = sum(rx)*dx
rx /= norm
norm = sum(ry)*dx
ry /= norm

lx, = ax.plot(x,rx*Lx, c='red')
ly, = ax.plot(y,ry, c='blue')
ax.set_xlim([min(min(x),min(y)), max(max(x),max(y))] )


axT = plt.axes([0.1,0.1,0.65,0.01])
slider = Slider(axT, "time", 0, Nfig-1, valinit=0, valstep=1)

axYval = plt.axes([0.1,0.15,0.65,0.01])
sliderYval = Slider(axYval, "Ymax", 0.1, 5, valinit=5, valstep=0.01)

def update(val):
    t = slider.val
    try:
        #r = np.loadtxt(dirname+"r{}.dat".format(int(t) ), delimiter=';')
        #ry = np.sum(r, axis = 0 )*dx
        #rx = np.sum(r, axis = 1 )*dy
        rx = np.loadtxt(dirname+"rx{}.dat".format(int(t) ), delimiter=';')
        ry = np.loadtxt(dirname+"ry{}.dat".format(int(t) ), delimiter=';')
        norm = sum(rx)*dx
        rx /= norm
        norm = sum(ry)*dy
        ry /= norm

        lx.set_ydata( rx*Lx )
        lx.set_linestyle('-')

        ly.set_ydata( ry )
        ly.set_linestyle('-')
    except:
        #r = np.loadtxt(dirname+"r.dat", delimiter=';')
        #ry = np.sum(r, axis = 0 )*dx
        #rx = np.sum(r, axis = 1 )*dy

        rx = np.loadtxt(dirname+"rx.dat", delimiter=';')
        ry = np.loadtxt(dirname+"ry.dat", delimiter=';')

        norm = sum(rx)*dx
        rx /= norm
        norm = sum(ry)*dy
        ry /= norm

        lx.set_ydata( rx*Lx )
        lx.set_linestyle(':')

        ly.set_ydata( ry )
        ly.set_linestyle(':')


    fig.canvas.draw_idle()

def updateYval(val):
	ymax = sliderYval.val

	ax.set_ylim([0.,ymax])

	fig.canvas.draw_idle()

slider.on_changed(update)
sliderYval.on_changed(updateYval)

if True:
    q =  1.0488
    kg = 1.
    v0 = 5.
    vp = 5.
    x0 = 0
    alpha = 40
    D = 1.
    
    def v(R, v0, vp, x0):
        return v0 + vp*np.sin(R*np.pi*2/Lx)

    R = np.linspace(-Lx/2, Lx/2, 1000)

    dR = R[1] - R[0]
    vlist = v(R,v0,vp,x0)

    rhoA = rhoR_analytical(vlist, dR,alpha,D,q)
    rho2A = rho2R_analytical(vlist, dR,alpha,D,q,kg)

    norm = sum(rhoA)*dR
    rhoA /= norm
    norm = sum(rho2A)*dR
    rho2A /= norm

    ax.plot(R,rhoA*Lx, label="Analytical")
    ax.plot(R,rho2A*Lx, label="Analytical2")

ax.legend()




plt.show()


