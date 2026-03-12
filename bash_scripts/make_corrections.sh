#!/bash

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_AN_test.root 2 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_AN_test.root 2 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_AN_test.root 2 &

#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.6_loose_mx.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../data/correctionFiles/corrections_10.2_AN_mid_mx.root 2 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.4_loose_mx.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_AN_mid_mx.root 2 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.2_loose_mx.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_AN_mid_mx.root 2 &


#./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.4/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.4_AN.root 2 &
#./corrections/calcCorrections /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root /volatile/clas12/users/jphelan/SIDIS/generator/clasdis/10.6/final_skims/gen_skim.root ../data/correctionFiles/corrections_10.6_AN.root 2 &


./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_3d_AN.root 3 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_3d_AN.root 3 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_3d_AN.root 3 &

#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_no_match_AN.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_no_match_AN.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_no_match_AN.root 0 &

#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_no_match_AN_4D.root 0 &

#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_pT.root 2 pT 4 0 1.6 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_pT.root 2 pT 4 0 1.6 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_pT.root 2 pT 4 0 1.6 &
#wait

#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_sector_pi.root 2 sector_pi 6 0 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_sector_pi.root 2 sector_pi 6 1 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_sector_pi.root 2 sector_pi 6 1 7 &

#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_sector_e.root 2 sector_e 6 0 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_sector_e.root 2 sector_e 6 0 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_sector_e.root 2 sector_e 6 0 7 &

#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_phi_q.root 2 phi_q 6 0 360 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_phi_q.root 2 phi_q 6 0 360 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_phi_q.root 2 phi_q 6 0 360 &


#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_comb_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_sector_e.root 2 sector_e 6 0 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_sector_e.root 2 sector_e 6 1 7 &
#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_sector_e.root 2 sector_e 6 1 7 &

wait

#python ../plotting/plotCorrections.py corrections_10.2_AN_test.root bin
#python ../plotting/plotCorrections.py corrections_10.2_AN_test.root acc

#python ../plotting/plotCorrections.py corrections_10.4_AN_test.root bin
#python ../plotting/plotCorrections.py corrections_10.4_AN_test.root acc

#python ../plotting/plotCorrections.py corrections_10.6_AN_test.root bin
#python ../plotting/plotCorrections.py corrections_10.6_AN_test.root acc

python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root acc

python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root acc

python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root bin
python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.4_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.4_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.6_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.6_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN_4D.root acc
#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN_4D.root bin