#/bash/!

for i in {1..1}
do
	clasdis --nmax 25000  --trig 25000 --beam 10.2 --z 0.275 --targ deuteron
	mv eventfiles/clasdisde.00.e10.200.emn0.75tmn.09.xs58.05nb.dis.0000.dat /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/lund/clasdis_$i.dat
done
