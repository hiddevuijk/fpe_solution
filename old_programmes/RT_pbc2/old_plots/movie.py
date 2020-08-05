import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

from scipy.integrate import simps
from sys import exit



dirname = "data/"
Nfig = 1000
frames = np.arange(0,Nfig,1)

fig, ax = plt.subplots()

axl = ax


x = np.loadtxt(dirname+"x.dat")
dx = x[1] - x[0]
y = np.loadtxt(dirname+"y.dat")
dy = y[1] - y[0]
xmin = x[0]
xmax = x[-1]
ymin = y[0]
ymax = y[-1]



time = np.loadtxt(dirname+"t.dat", delimiter=';')
dt = time[2] - time[1]

rx = np.loadtxt(dirname+"rx.dat", delimiter=';')
rxmax = 1.1*max(rx)

rx = np.loadtxt(dirname+"rx0.dat", delimiter=';')
ry = np.loadtxt(dirname+"ry0.dat", delimiter=';')

lx, = axl.plot(x,rx, c='red')
#ly, = ax.plot(y,ry, c='blue')
def init_func():
	axl.clear()
	axl.set_xlabel("x")
	axl.set_ylabel(r"$\rho(x)$")
	axl.set_ylim([0,rxmax])
	axl.set_title("t={:.2f}".format(time[0]) )


def update_plot(i):
	t = i
	
	try:
		axl.clear()
		rx = np.loadtxt(dirname+"rx{}.dat".format(int(t) ), delimiter=';')
		#ry = np.loadtxt(dirname+"ry{}.dat".format(int(t) ), delimiter=';')
		axl.plot(x,rx,c='red')
		axl.set_ylim([0,rxmax])
		axl.set_title("t={:.2f}".format(t*dt) )
		axl.set_xlabel("x")
		axl.set_ylabel(r"$\rho(x)$")
		axl.set_ylim([0,rxmax])

		

	except:
		rx = np.loadtxt(dirname+"rx.dat", delimiter=';')
		#ry = np.loadtxt(dirname+"ry.dat", delimiter=';')

		axl.plot(x,rx,c='red')
		axl.set_ylim([0,rxmax])
		axl.set_title("t={:.2f} steady state".format(t*dt) )
		axl.set_xlabel("x")
		axl.set_ylabel(r"$\rho(x)$")
		axl.set_ylim([0,rxmax])





anim = FuncAnimation(fig, update_plot, frames=frames, init_func=init_func, interval=20)

#plt.show()
anim.save("anim.mp4", dpi=150, fps=100, writer='ffmpeg')

