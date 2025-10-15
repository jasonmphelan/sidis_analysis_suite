#!/bin/bash

for dir in /mss/clas12/rg-b/production/recon/fall2019/torus+1/pass2/v1/dst/recon/*/;
do
	jcache get $dir*
done
