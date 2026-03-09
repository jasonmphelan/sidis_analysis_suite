#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
#import numpy as np
import time

Ebeam = 10.6
# Create a list with the ranges:

for x in range(0,100):
	n = x
	#n = runList[x]
	time.sleep(0.5)			
	command="""#!/bin/sh 
#SBATCH --job-name=rotationSym_{0}_{1}
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --mem-per-cpu=500
#SBATCH -t60
#SBATCH --constraint=el9
#SBATCH --error=/farm_out/jphelan/rotation/err_rot_{0}_{1}.txt
#SBATCH --output=/volatile/clas12/users/jphelan/SIDIS/cut_sensitivity/farm_out/out_{0}_{1}.txt
time ./../rhoAnalysis/rotateRhoSym /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_{1}.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_{1} {1} 1 0 acceptanceMap_allE_final.root {0} 
    """.format( n, Ebeam)	
	command = command.replace('\t', '')
	print(command)
	p=Popen(args=["sbatch"],stdin=PIPE)
	p.communicate(command.encode())