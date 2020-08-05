import numpy as np
from sys import exit
from os import mkdir
from shutil import copy

# box size 
Lx = 10.
Nx = 100
Ny = 100


# potential
kList = [0., 4., 8.]

#swims tumble
v0 = 40.
vp = 8.
vx0 = 0
alpha = 20.

#diffusion
gamma = 1.
GammaList = [2/5., 3/5., 4/5., 1., 5./4, 5./3, 5./2, 5.]
temp = 1.

#integration
time = 0.1
nc = 0.00075
saveXY = 0

Nprint = 100
Nsave = 1000

epsilon = 0.0001
Nss = 10000


program = 'go.exe'
dirnames = open('dirnames.txt','w')

for ki, k in enumerate(kList):
	if k == 0: Ly = 1.
	else:      Ly = 5.
	dy = Ly/Ny

	for Gi, G in enumerate(GammaList):
			dt = nc * G*dy*dy/temp
			print(dt)

			dirname = "k{}_G{}".format(ki, Gi)
			dirnames.write(dirname+'\n')
			mkdir(dirname)

			copy(program,dirname+'/'+program)
			infile = open(dirname+"/input.txt",'w')

			infile.write("# boxsize \n")
			infile.write("Lx = {:.10f}\n".format(Lx) )
			infile.write("Nx = {}\n".format(Nx) )
			infile.write("Ly = {:.10f}\n".format(Ly) )
			infile.write("Ny = {}\n".format(Ny) )

			infile.write("\n# potential \n")
			infile.write("k = {:.5f}\n".format(k) )

		
			infile.write("\n# swims tumble \n")
			infile.write("v0 = {:.10f}\n".format(v0) )
			infile.write("vp = {:.10f}\n".format(vp) )
			infile.write("vx0 = {:.10f}\n".format(vx0) )
			infile.write("alpha = {:.20f}\n".format(alpha) )


			infile.write("\n# diffusion\n")
			infile.write("gamma = {:.10f}\n".format(gamma) )
			infile.write("Gamma = {:.10f}\n".format(G) )
			infile.write("temp = {:.10f}\n".format(temp) )


			infile.write("\n# integration\n")
			infile.write("dt = {:.16f}\n".format(dt) )
			infile.write("time = {:.5f}\n".format(time) )
			infile.write("Nprint = {}\n".format(Nprint) )
			infile.write("Nsave = {}\n".format(Nsave) )
			infile.write("saveXY = {}\n".format(saveXY) )

			infile.write("epsilon = {:.10f}\n".format(epsilon) )
			infile.write("Nss = {}\n".format(Nss) )
			
		
			
			infile.close()


			runfile = open(dirname+"/run.sh",'w')
			runfile.write("#!/bin/bash\n")
			runfile.write("#SBATCH -J " + dirname+'\n')
			runfile.write("#SBATCH -n 1"+'\n')
			runfile.write("#SBATCH -e out.%J"+'\n')
			runfile.write("#SBATCH -o out.%j"+'\n')
			runfile.write("#SBATCH -p all"+'\n')
			runfile.write("time ./"+program+'\n')
			runfile.close()

dirnames.close()
