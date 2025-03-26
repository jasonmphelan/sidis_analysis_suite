#!/bin/csh

echo SCRIPTHEADER START:  `date +%s`
printf 'Job is running on node: '; /bin/hostname
printf 'Job submitted by: jphelan'
echo Running directory: `pwd`
set submissionID=$1
set sjobID=$2
echo Directory `pwd` content before starting submissionID $submissionID':'
ls -l
if ($? != 0) then
  echo ls failure
  exit 211
endif

# Run Script Header
# -----------------
set echo
# limit core size
# ulimit -c 10

# clearing module environment, as suggested by OSG, #67051
unsetenv ENABLE_LMOD
unsetenv _LMFILES_
unsetenv LMOD_ANCIENT_TIME
unsetenv LMOD_arch
unsetenv LMOD_CMD
unsetenv LMOD_COLORIZE
unsetenv LMOD_DEFAULT_MODULEPATH
unsetenv LMOD_DIR
unsetenv LMOD_FULL_SETTARG_SUPPORT
unsetenv LMOD_PACKAGE_PATH
unsetenv LMOD_PKG
unsetenv LMOD_PREPEND_BLOCK
unsetenv LMOD_SETTARG_CMD
unsetenv LMOD_SETTARG_FULL_SUPPORT
unsetenv LMOD_sys
unsetenv LMOD_SYSTEM_DEFAULT_MODULES
unsetenv LMOD_VERSION
unsetenv LOADEDMODULES
unsetenv MODULEPATH
unsetenv MODULEPATH_ROOT
unsetenv MODULESHOME

# Exit if cvmfs file not found, source it if found
# ------------------------------------------------

source /etc/profile.d/modules.csh
set cvmfsPath = /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/
set cvmfsSetupFile = $cvmfsPath/setup.csh
if (-f $cvmfsSetupFile ) then
		echo $cvmfsSetupFile exists, sourcing it with path $cvmfsPath
		source $cvmfsSetupFile $cvmfsPath
		set cvmfsSetupFile2 = /cvmfs/oasis.opensciencegrid.org/jlab/geant4/ceInstall/geant4_cvmfs.csh 
		if (-f $cvmfsSetupFile2 ) then
			echo $cvmfsSetupFile2 exists, sourcing it
			source $cvmfsSetupFile2 
		else
			echo CVMFS ERROR $cvmfsSetupFile2 does not exist. Exiting
			exit 202
		endif
endif
else
		echo CVMFS ERROR $cvmfsSetupFile does not exist. Exiting
		exit 202
endif
if ( ! -f /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/gemc/5.10/rgb_spring2019.gcard ) then
	echo gcard not found, exiting
	exit 241
endif
if ( ! -f /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/coatjava/10.0.7/rgb_spring2019.yaml ) then
	echo yaml not found, exiting
	exit 242
endif

module unload gemc
module unload coatjava
module unload jdk
module unload root
module unload mcgen
module load gemc/5.10
module load sqlite/5.10

# TODO: RCDB_CONNECTION currently not used. When fixed, remove this line.
setenv RCDB_CONNECTION mysql://null

module avail
module load coatjava/10.0.7
module load jdk/17.0.2
#module load root/6.28.06
module load mcgen/3.12

echo ROOT Version: 6.28.06
echo MCGEN Version: 3.12
echo SQLITE Version: 5.10
echo GEMC Version: 5.10
echo JDK Version: 17.0.2
echo COATJAVA Version: 10.0.7




echo
echo
echo SCRIPTHEADER END:  `date +%s`
echo

# End of Run Script Header
# ------------------------


# Generator
# ---------

# saving date for bookmarking purposes:
echo GENERATOR START:  `date +%s`

# generate-seeds.py is in the path
generate-seeds.py generate


set seed = `generate-seeds.py read --row 1`
echo Generator seed from generate-seeds, row 1: $seed

echo
echo Running 10000 events with generator clasdis with options: --beam 10.6 --targ deuteron --z 0.275 
echo Generator:
which clasdis
echo
clasdis --trig 10000 --docker --beam 10.6 --targ deuteron --z 0.275 --seed $seed

if ($? != 0) then
  echo GENERATOR ERROR '>'clasdis'<' failed.
  exit 203
endif

# removing any root file that was created
rm -f *.root

echo
printf 'Events Generator Completed on: '; /bin/date
echo
echo 'Directory Content After Generator:'
ls -l
if ($? != 0) then
  echo ls failure
  exit 211
endif
echo

echo GENERATOR END:  `date +%s`

# End of Run Generator
# ---------------------


# Run GEMC
# --------

echo GEMC START:  `date +%s`

echo
echo GEMC executable: 
which gemc
echo gcard: /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/gemc/5.10/rgb_spring2019.gcard
echo
ls ../

gemc -USE_GUI=0 -N=10000 /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/gemc/5.10/rgb_spring2019.gcard  -INPUT_GEN_FILE='lund, clasdis.dat'   -RANDOMIZE_LUND_VZ='-3.0*cm, 2.5*cm, reset '  -BEAM_SPOT='0.0*mm, 0.0*mm, 0.0*mm, 0.0*mm, 0*deg, reset '          -RASTER_VERTEX='0.0*cm, 0.0*cm, reset '       -SCALE_FIELD='binary_torus,    -1.00'   -SCALE_FIELD='binary_solenoid, -1.00'   -OUTPUT='hipo, gemc.hipo'   -INTEGRATEDRAW='*'   | sed '/G4Exception-START/,/G4Exception-END/d'  
if ($? != 0) then
	echo gemc failed
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 204
endif
# removing generated events file
rm -f *.dat

echo
printf 'GEMC Completed on: '; /bin/date
echo
echo 'Removing LUND file'
rm -f clasdis.dat
echo
echo 'Directory Content After GEMC:'
ls -l
if ($? != 0) then
	echo ls failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 211
endif
echo

echo GEMC END:  `date +%s`

# End of GEMC
# -----------

echo Gemc 5.1 or later has hipo output, no need to run evio2hipo
  
# Run de-noiser
# -------------

echo DE-NOISING START:  `date +%s`

$DRIFTCHAMBERS/install/bin/denoise2.exe  -i gemc.hipo -o gemc_denoised.hipo -t 1 -n $DRIFTCHAMBERS/denoising/code/network/cnn_autoenc_0f_112.json 

if ($? != 0) then
	echo de-noiser failed.
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 230
endif

echo Removing original hipo file gemc.hipo
rm -f gemc.hipo
echo 'Directory Content After de-noiser:'
ls -l
if ($? != 0) then
	echo ls failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 211
endif

echo DE-NOISING END:  `date +%s`

# End of de-noiser
# ----------------




# Run Reconstruction
# ------------------

echo RECONSTRUCTION START:  `date +%s`

echo content of yaml file /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/coatjava/10.0.7/rgb_spring2019.yaml:
cat /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/coatjava/10.0.7/rgb_spring2019.yaml

echo
df /cvmfs/oasis.opensciencegrid.org && df . && df /tmp
if ($? != 0) then
	echo df failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 213
endif

echo
echo executing: recon-util -y /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/coatjava/10.0.7/rgb_spring2019.yaml -i gemc_denoised.hipo -o recon.hipo
recon-util -y /cvmfs/oasis.opensciencegrid.org/jlab/hallb/clas12/sw/noarch/clas12-config/prod/coatjava/10.0.7/rgb_spring2019.yaml -i gemc_denoised.hipo -o recon.hipo
if ($? != 0) then
	echo recon-util failed.
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 207
endif
df /cvmfs/oasis.opensciencegrid.org && df . && df /tmp
if ($? != 0) then
	echo df failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 213
endif

echo
printf 'recon-util Completed on: '; /bin/date
echo
echo 'Directory Content After recon-util:'
ls -l
if ($? != 0) then
	echo ls failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 211
endif

echo executing: hipo-utils -test recon.hipo
hipo-utils -test recon.hipo
if ($? != 0) then
	echo hipo-utils failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 214
endif

if (`stat -L -c%s recon.hipo` < 100) then
	echo hipo size failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 215
endif

echo
echo RECONSTRUCTION END:  `date +%s`
echo

# End of Reconstruction
# ---------------------


# Removing Unnecessary Files and Creating DST if selected
# -------------------------------------------------------


echo
echo Creating the DST
echo
hipo-utils -filter -b 'RUN::*,RAW::epics,RAW::scaler,HEL::flip,HEL::online,REC::*,RECFT::*,MC::RecMatch,MC::GenMatch,MC::Particle,MC::User,MC::Header,MC::Lund,MC::Event' -merge -o dst.hipo recon.hipo
set outputFileName=''$submissionID'-'$sjobID'.hipo'
echo submissionID is $submissionID
echo sjobID is $sjobID
echo outputFileName is $outputFileName
echo Moving the DST to the output file $outputFileName
mv dst.hipo $outputFileName
if ($? != 0) then
  echo hipo-utils failed, removing data files and exiting
  rm -f *.hipo *.evio
  exit 208
endif

echo
printf 'DST Completed on: '; /bin/date
echo
echo 'Directory Content After DST:'
ls -l
if ($? != 0) then
	echo ls failure
	echo removing data files and exiting
	rm -f *.hipo *.evio
	exit 211
endif
echo
echo
echo Simulation + Reconstruction Successfully Completed on: `date +%s` 


echo Additional cleanup
rm -f core* *.gcard
rm -f recon.hipo gemc.hipo gemc.merged.hipo gemc_denoised.hipo 
rm -f run.sh nodeScript.sh condor_exec.exe
rm -f RNDMSTATUS random-seeds.txt clasdis.dat
rm -f gemc.evio

echo
echo nodeScript run completed.
echo 'Final Directory Content:'
ls -l
echo

# End of Run Script Footer
# ------------------------

	
