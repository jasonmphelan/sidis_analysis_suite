#!/bash

#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10820.root ../data/correctionFiles/corrections_10.2_10820.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10637.root ../trees/final_skims/GEMC/gen_skim_10637.root ../data/correctionFiles/corrections_10.2_10637.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10892.root ../data/correctionFiles/corrections_10.4_10892.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../data/correctionFiles/corrections_10.2.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6.root 1 0 0 matchCut2D_10.root 0 &

#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_no_map.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

./../build/corrections/calcCorrections ../data/correctionFiles/corrections_allE.root 1 0 5 5 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &


#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_3d_AN.root 3 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_3d_AN.root 3 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_3d_AN.root 3 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_rga_pT.root 1 1 pT 4 0 1.2 1 1 ../trees/final_skims/GEMC/final_skim_10927_rga.root ../trees/final_skims/GEMC/gen_skim_10927_rga.root matchCut2D_10.root &
./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_pT.root 1 0 pT 4 0 1.2 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.4_pT.root 1 0 pT 4 0 1.2 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.6_pT.root 1 0 pT 4 0 1.2 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_AN.root 1 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_AN.root 1 0 0 matchCut2D_10.root 0 &

#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10.2_merged.root ../data/correctionFiles/corrections_10.2_3d_AN.root 3 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_3d_AN.root 3 0 0 matchCut2D_10.root 0 &
#./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_3d_AN.root 3 0 0 matchCut2D_10.root 0 &


#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_no_match.root 0 0 pT 4 0 1.6 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &

#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_pT.root 1 0 pT 4 0 1.6 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.4_pT.root 1 0 pT 4 0 1.6 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.6_pT.root 1 0 pT 4 0 1.6 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_Mx.root 1 0 Mx 4 1.7 2.9 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.4_Mx.root 1 0 Mx 4 1.7 2.9 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.6_Mx.root 1 0 Mx 4 1.7 2.9 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &



#./../build/corrections/calcCorrections4D ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10.2_merged.root ../data/correctionFiles/corrections_10.2_pT.root 1 0 pT 4 0 1.6 &
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

#python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.2_3d_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.4_3d_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.6_3d_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.4_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.4_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.6_no_match_AN.root bin
#python ../plotting/plotCorrections.py corrections_10.6_no_match_AN.root acc

#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN_4D.root acc
#python ../plotting/plotCorrections.py corrections_10.2_no_match_AN_4D.root bin