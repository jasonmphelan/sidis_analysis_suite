#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
import numpy as np
import time
import sys


Eindex = int(sys.argv[1])
Npindex = int(sys.argv[2])

E_options = [10.2, 10.4, 10.6]
E = E_options[Eindex]
Ebeam = E_options[Eindex]

Np_options = [50, 150, 250, 350, 450]
Np = Np_options[Npindex]

xsv = 'tagged_full'
xpv = 'PRC'
variation='';
suffix='_lowAlphaS'

#exe='/work/clas12/users/tkutz/gemc/clas12Tags/4.4.1/source/gemc'
#gcard='/work/clas12/users/tkutz/gemc/clas12Tags/gcards/band-clas12{0}.gcard'.format(variation)

exe='gemc'
gcard='/work/clas12/users/tkutz/gemc/clas12-config/gemc/5.10/rgb_spring2019.gcard'
lunddir='/volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/neutron_'

#reruns = np.loadtxt("rerun_{0}MeV{1}.dat".format(Np,suffix), dtype=int)
#for n in reruns:

#for n in np.arange(11,1001):

for n in [442, 670, 908]:
				
	time.sleep(0.5)			
	command="""#!/bin/bash -l
#SBATCH --job-name={0}_{1}_{8}GeV_{10}_{7:05d}
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=1250
#SBATCH -t720
#SBATCH --constraint=el9
#SBATCH --error=/farm_out/tkutz/gemc/err_gemc_{0}_{1}{6}{9}_{8}GeV_{10}MeV_{7:05d}.txt
#SBATCH --output=/volatile/clas12/users/tkutz/GEMC/farm_out/gemc/tagged_full/{8}/out_gemc_{0}_{1}{6}{9}_{8}GeV_{10}MeV_{7:05d}.txt
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12
module switch gemc/5.10
module switch sqlite/5.10
export GEMC_DATA_DIR="/work/clas12/users/tkutz/gemc/clas12Tags/clas12Tags-5.10"
export CCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/ccdb_05-12-2024.sqlite"
export RCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/rcdb_2024-06-18.sqlite"
time {2} {3} -USE_GUI=0 -N=25000 -INPUT_GEN_FILE="LUND, {4}/lund_{0}_{1}_E{8}GeV_Np{10}MeV{6}_{7:05d}.dat" -OUTPUT="hipo, {5}/gemc_{0}_{1}_E_{8}GeV{9}_Np_{10}MeV{6}_{7:05d}.hipo"
	""".format(xsv, xpv, exe, gcard, lunddir, outdir, suffix, n, Ebeam, variation, Np)	
	command = command.replace('\t', '')
	print(command)
	p=Popen(args=["sbatch"],stdin=PIPE)
	p.communicate(command.encode())
