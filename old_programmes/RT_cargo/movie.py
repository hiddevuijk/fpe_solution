import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from scipy.integrate import simps
from sys import exit

dirname = "data/"
Nfig = 5000

fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.2)



x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]

r = np.loadtxt(dirname+"r0.dat", delimiter=';')
ry = np.sum(r, axis = 0 )*dx
rx = np.sum(r, axis = 1 )*dy

#normy = simps(ry,y)
#ry /= normy
#
#normx = simps(rx,x)
#rx /= normx

lx, = ax.plot(x,rx, c='red')
ly, = ax.plot(y,ry, c='blue')
ryavg = simps(y*ry, y)
lavgy = ax.axvline(ryavg, c='blue')
rxavg = simps(x*rx, x)
lavgx = ax.axvline(rxavg, c='red')
ax.set_xlim([min(min(x),min(y)), max(max(x),max(y))] )


axT = plt.axes([0.1,0.1,0.65,0.01])
slider = Slider(axT, "time", 0, Nfig-1, valinit=0, valstep=1)

axYval = plt.axes([0.1,0.15,0.65,0.01])
sliderYval = Slider(axYval, "Ymax", 0.1, 5, valinit=5, valstep=0.01)

def update(val):
	t = slider.val
	r = np.loadtxt(dirname+"r{}.dat".format(int(t) ), delimiter=';')
	ry = np.sum(r, axis = 0 )*dx
	rx = np.sum(r, axis = 1 )*dy



	lx.set_ydata( rx )
	ly.set_ydata( ry )

	ryavg = simps(y*ry, y)
	lavgy.set_xdata(ryavg)

	rxavg = simps(x*rx, x)
	lavgx.set_xdata(rxavg)

	fig.canvas.draw_idle()

def updateYval(val):
	ymax = sliderYval.val

	ax.set_ylim([0.,ymax])

	fig.canvas.draw_idle()

slider.on_changed(update)
sliderYval.on_changed(updateYval)

plt.show()


