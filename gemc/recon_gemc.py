#!/usr/bin/env python3

from subprocess import Popen, PIPE
import subprocess
import numpy as np
import sys
import time 

Eindex = int(sys.argv[1])

E_options = [10.2, 10.4, 10.6]
E = E_options[Eindex]
Ebeam = E_options[Eindex]


xsv = 'tagged_full'
xpv = 'PRC'
variation='';
suffix='_lowAlphaS'
nuc = 'proton'

decdir='/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/{0}/hipo'.format(Ebeam)
recdir='/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/{0}/reco_out'.format(Ebeam)
yamlfile='/work/clas12/users/tkutz/gemc/clas12-config/coatjava/10.0.2/rgb_spring2019.yaml'

#reruns = np.loadtxt('rerecon_{0}MeV{1}.dat'.format(Np,suffix),dtype=int)
#for n in reruns:

for n in range(5000):

    time.sleep(0.25)			

    decfile = '{0}/{1}_{2}.hipo'.format(decdir,nuc,n+1)
    recfile = '{0}/{1}_{2}.hipo'.format(recdir,nuc,n+1)
    errfile = '/farm_out/jphelan/clasdis/err_gemc_{0}_{1}.txt'.format(nuc, n+1)
    outfile = '/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/farm_out/out_gemc_{0}_{1}.txt'.format(nuc, n+1)
	
    command="""#!/bin/bash -l
#SBATCH --job-name=cook_gemc
#SBATCH --account=clas12
#SBATCH -p production
#SBATCH --constraint=el9
#SBATCH --export=NONE
#SBATCH --mem-per-cpu=1500
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
