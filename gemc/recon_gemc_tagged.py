#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
import numpy as np
import sys
import time 

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

decdir='/volatile/clas12/users/tkutz/GEMC/decode/{2}/{0}/{1}'.format(xsv,xpv,Ebeam)
recdir='/volatile/clas12/users/tkutz/GEMC/recon/{2}/{0}/{1}'.format(xsv,xpv,Ebeam)
yamlfile='/work/clas12/users/tkutz/gemc/clas12-config/coatjava/10.0.2/rgb_spring2019.yaml'

#reruns = np.loadtxt('rerecon_{0}MeV{1}.dat'.format(Np,suffix),dtype=int)
#for n in reruns:

for n in np.arange(11,1001):

	time.sleep(0.25)			

	decfile = '{0}/gemc_{1}_{2}_E_{3}GeV_Np_{5}MeV{6}_{4:05d}.hipo'.format(decdir, xsv, xpv, Ebeam, n, Np, suffix)
	recfile = '{0}/recon_{1}_{2}_E_{3}GeV_Np_{5}MeV{6}_{4:05d}.hipo'.format(recdir, xsv, xpv, Ebeam, n, Np, suffix)
	errfile = '/farm_out/tkutz/cook_gemc/{0}/err_cook_{0}_{1}_E_{2}GeV_Np_{4}MeV{5}_{3:05d}.txt'.format(xsv, xpv, Ebeam, n, Np, suffix) 
	outfile = '/farm_out/tkutz/cook_gemc/{0}/out_cook_{0}_{1}_E_{2}GeV_Np_{4}MeV{5}_{3:05d}.txt'.format(xsv, xpv, Ebeam, n, Np, suffix) 
	
	command="""#!/bin/bash -l
#SBATCH --job-name=cook_gemc
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --constraint=el9
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=2000
#SBATCH -t600
#SBATCH --error=%s
#SBATCH --output=%s
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12 tmpfs sqlite/5.9
module switch coatjava/10.0.2
module switch ccdb/1.99.1
module switch rcdb/1.99.0
export CCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/ccdb_05-12-2024.sqlite"
export RCDB_CONNECTION="sqlite:////work/clas12/users/tkutz/gemc/database/rcdb_2024-06-18.sqlite"
recon-util -i %s -o %s -y %s 	
	""" % (errfile, outfile, decfile, recfile, yamlfile)
	command = command.replace('\t', '')
	print(command)
	p=Popen(args=["sbatch"],stdin=PIPE)
	p.communicate(command.encode())
