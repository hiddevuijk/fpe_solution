#!/bin/bash

while read dirname
do
	cd /beetmp/hidde/density/rho01/$dirname
	sbatch run.sh
done < dirnames.txt



