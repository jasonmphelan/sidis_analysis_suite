#!/bin/bash

declare -i var1=0

for dir in /cache/clas12/rg-b/production/recon/fall2019/torus+1/pass2/v1/dst/recon/*/;
do
	echo $var1
	echo /volatile/clas12/users/jphelan/SIDIS/background/charge_symm_$var1.hipo
	for sec in {7,8,9,10,11,12}
	do
		trigger-filter -b $sec -o /volatile/clas12/users/jphelan/SIDIS/data/background_hipo/run_10.4_out_${sec}_${var1}.hipo $dir*
	done
	#jcache get $dir*
	((var1++))
done