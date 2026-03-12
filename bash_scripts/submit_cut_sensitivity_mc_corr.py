#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
#import numpy as np
import time

Ebeam = 10.2

for x in range(0, 300):
	n = x + 1 
	#n = runList[x]
	time.sleep(0.5)			
	command="""#!/bin/sh 
#SBATCH --job-name=cut_sensitivity_{0}_{1}
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --mem-per-cpu=250
#SBATCH -t600
#SBATCH --constraint=el9
#SBATCH --error=/farm_out/jphelan/cut_sensitivity/err_cs_{0}_{1}.txt
#SBATCH --output=/volatile/clas12/users/jphelan/SIDIS/cut_sensitivity/farm_out/out_{0}_{1}.txt
time ./cutAnalysis/cutSensitivity 10 0 0 0 0 /volatile/clas12/users/jphelan/SIDIS/cut_sensitivity/ratio_no_vz_{1}_{0}
    """.format( n, Ebeam)	
	command = command.replace('\t', '')
	print(command)
	p=Popen(args=["sbatch"],stdin=PIPE)
	p.communicate(command.encode())
#SBATCH --export=NONE
#module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
#module load clas12
#export CCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/ccdb_05-12-2024.sqlite"
#export RCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/rcdb_2024-06-18.sqlite"
#time {2} {3} -USE_GUI=0 -N=25000 -INPUT_GEN_FILE="LUND, {4}{7}.dat" -RANDOMIZE_LUND_VZ='-3.0*cm, 2.5*cm, reset '  -BEAM_SPOT='0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0*deg, reset '          -RASTER_VERTEX='0.0*cm, 0.0*cm, reset '       -SCALE_FIELD='binary_torus,    -1.00'   -SCALE_FIELD='binary_solenoid, -1.00' -OUTPUT="hipo, /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/{8}/hipo/proton_{7}.hipo"
#export GEMC_DATA_DIR="/work/clas12/users/jphelan/GEMC_DATA/clas12Tags-5.10"
