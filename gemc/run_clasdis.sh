#/bash/!

for i in {1..2000}
do
	clasdis --nmax 5000  --trig 5000 --beam 10.2 --z 0.275 --targ deuteron
	#clasdis --nmax 5000 --trig 5000 --beam 10.2 --z 0 --zpos -3 --targ proton
	#clasdis --nmax 5000 --trig 5000 --beam 10.2 --z 0 --zpos -3 --targ neutron
	mv eventfiles/clasdisde.00.e10.200.emn0.75tmn.09.xs58.05nb.dis.0000.dat  /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/deuteron_$i.dat
	#mv eventfiles/clasdispr.00.e10.200.emn0.75tmn.09.xs65.59nb.dis.0000.dat /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/proton_$i.dat
	#mv eventfiles/clasdisne.00.e10.200.emn0.75tmn.09.xs50.51nb.dis.0000.dat /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/neutron_$i.dat
done
