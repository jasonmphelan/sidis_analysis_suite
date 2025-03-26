#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
#import numpy as np
import time

Ebeam = 10.2
xsv = 'deuteron'
xpv = 'dis'
suffix='';

#Note, everything works with version 4.4.1
gcard='/cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/gemc/5.10/rgb_spring2019.gcard'
exe='gemc'
gcard='/work/clas12/users/jphelan/GEMC_DATA/rgb_spring2019.gcard'
lunddir='/volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/clasdis_'


for x in range(1):
	n = x + 1 
	#n = runList[x]
	time.sleep(0.5)			
	command="""#!/bin/sh 
#SBATCH --job-name={0}_{1}_{8}GeV_{7:05d}
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --mem-per-cpu=1500
#SBATCH -t720
#SBATCH --constraint=el9
#SBATCH --error=/farm_out/jphelan/clasdis/err_gemc_{0}_{7}.txt
#SBATCH --output=/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/farm_out/out_gemc_{0}_{7}.txt
time {2} {3} -USE_GUI=0 -N=25000 -INPUT_GEN_FILE="LUND, {4}{7}.dat" -RANDOMIZE_LUND_VZ='-3.0*cm, 2.5*cm, reset '  -BEAM_SPOT='0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0*deg, reset '          -RASTER_VERTEX='0.0*cm, 0.0*cm, reset '       -SCALE_FIELD='binary_torus,    -1.00'   -SCALE_FIELD='binary_solenoid, -1.00' -OUTPUT="hipo, /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/{8}/hipo/clasdis_{7}.hipo"
	""".format(xsv, xpv, exe, gcard, lunddir, 0, suffix, n, Ebeam)	
	command = command.replace('\t', '')
	print(command)
	p=Popen(args=["sbatch"],stdin=PIPE)
	p.communicate(command.encode())
#SBATCH --export=NONE
#module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
#module load clas12
#export CCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/ccdb_05-12-2024.sqlite"
#export RCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/rcdb_2024-06-18.sqlite"
