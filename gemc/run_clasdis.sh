#/bash/!

for i in {1..5000}
do
	#clasdis --nmax 5000  --trig 5000 --beam 10.2 --z 0 --zpos -3 --targ deuteron --q 1.75 10
	clasdis --nmax 5000 --trig 5000 --beam 10.2 --z 0 --zpos -3 --targ proton --q 1.75 10
	clasdis --nmax 5000 --trig 5000 --beam 10.2 --z 0 --zpos -3 --targ neutron --q 1.75 10
	#mv eventfiles/clasdisde.00.e10.200.emn0.75tmn.09.xs58.05nb.dis.0000.dat  /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/deuteron_$i.dat
	mv eventfiles/clasdispr.00.e10.200.emn0.75tmn.09.xs22.32nb.dis.0000.dat /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/proton_$i.dat
	mv eventfiles/clasdisne.00.e10.200.emn0.75tmn.09.xs15.34nb.dis.0000.dat /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/neutron_$i.dat
done
