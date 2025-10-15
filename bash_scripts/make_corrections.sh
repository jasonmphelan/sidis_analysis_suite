#!/bash

./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.2_AN.root 2 &
./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.4/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.4_AN.root 2 &
./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.6/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.6_AN.root 2 &

./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.2/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.2_3d_AN.root 3 &
./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.4/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.4_3d_AN.root 3 &
./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.6/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.6_3d_AN.root 3 &

wait

python ../plotting/plotCorrections.py corrections_10.2_AN.root bin
python ../plotting/plotCorrections.py corrections_10.2_AN.root acc

python ../plotting/plotCorrections.py corrections_10.4_AN.root bin
python ../plotting/plotCorrections.py corrections_10.4_AN.root acc

python ../plotting/plotCorrections.py corrections_10.6_AN.root bin
python ../plotting/plotCorrections.py corrections_10.6_AN.root acc

python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root acc

python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root acc

python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root acc
