#!/bash

################
# Corrections
################

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_3d.root 3 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_3d.root 3 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_3d.root 3 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_pT.root 1 0 pT 4 0 1.2 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.4_pT.root 1 0 pT 4 0 1.2 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.6_pT.root 1 0 pT 4 0 1.2 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_pT_no_match.root 0 0 pT 5 0 1.5 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.4_pT_no_match.root 0 0 pT 5 0 1.5 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_10.root &
#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.6_pT_no_match.root 0 0 pT 5 0 1.5 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#./../build/corrections/calcCorrections4D ../data/correctionFiles/corrections_10.2_rga_pT.root 1 1 pT 4 0 1.2 1 1 ../trees/final_skims/GEMC/final_skim_10927_rga.root ../trees/final_skims/GEMC/gen_skim_10927_rga.root matchCut2D_10.root &


wait


#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.2_pT.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.2_pT_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.4_pT.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.4_pT_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.6_pT.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.6_pT_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" & 

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.2_pT_no_match.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.2_pT_no_match_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.4_pT_no_match.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.4_pT_no_match_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.6_pT_no_match.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.6_pT_no_match_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" & 

papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE_4D ../data/correctionFiles/corrections_10.2_pT_no_match.root -p OUTPUT_FILE_4D ../data/correctionFiles/corrections_10.2_no_match_pT_fit.root -y "OUTPUT_FILE: null " -y "FILE: ../data/correctionFiles/corrections_10.2_AN.root" & 

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_3d.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_3d_fit.root -y "OUTPUT_FILE_4D: null " -y "FILE_4D: ../data/correctionFiles/corrections_10.2_pT.root" &


wait

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_1.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_1.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_1.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_1.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_1.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_1.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_5.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_5.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_5.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_5.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_5.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_10.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_15.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_15.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_15.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_15.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_15.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_15.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_20.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_20.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_20.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_20.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_20.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_20.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_25.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_25.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_25.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_25.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_25.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_25.root &

#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.2_AN_50.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10820.root ../trees/final_skims/GEMC/gen_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_10820.root matchCut2D_50.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.4_AN_50.root 1 0 2 2 ../trees/final_skims/GEMC/final_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10892.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10892.root matchCut2D_50.root &
#../build/corrections/calcCorrections ../data/correctionFiles/corrections_10.6_AN_50.root 1 0 1 1 ../trees/final_skims/GEMC/final_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root matchCut2D_50.root &

wait 


#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_1.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_1_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_1.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_1_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_1.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_1_fit.root &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_5.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_5_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_5.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_5_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_5.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_5_fit.root &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_15.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_15_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_15.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_15_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_15.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_15_fit.root &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_20.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_20_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_20.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_20_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_20.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_20_fit.root &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_25.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_25_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_25.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_25_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_25.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_25_fit.root &

#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.2_AN_50.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.2_AN_50_fit.root &
#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.4_AN_50.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.4_AN_50_fit.root &#papermill ../plotting/plotCorrections3D.ipynb /dev/null -p FILE ../data/correctionFiles/corrections_10.6_AN_50.root -p OUTPUT_FILE ../data/correctionFiles/corrections_10.6_AN_50_fit.root &

wait

#################
#ANALYSIS RATIOS
#################

#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_mc.root 1 3 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 0 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 0 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d.root 1 0 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 0 1 0 matchCut2D_10.root & 

#################
#Energy RATIOS
#################
#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/final_note/ratios_2d_rho_10.2.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 10.4 ../histograms/final_note/ratios_2d_rho_10.4.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 10.6 ../histograms/final_note/ratios_2d_rho_10.6.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root & 

wait

#################
#Acc Match RATIOS
#################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_1.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_1_fit.root 1.25 1 0 matchCut2D_1.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_5.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_5_fit.root 1.25 1 0 matchCut2D_5.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_15.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_15_fit.root 1.25 1 0 matchCut2D_15.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_20.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_20_fit.root 1.25 1 0 matchCut2D_20.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_25.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_25_fit.root 1.25 1 0 matchCut2D_25.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_k_50.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_50_fit.root 1.25 1 0 matchCut2D_50.root & 

wait
#################
#3D RATIOS
#################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_3d.root 3 0 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_3d_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_3d_k.root 3 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_3d_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_3d_k_no_mc.root 3 4 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_3d_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_3d_k_aligned.root 4 4 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_3d_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_qInt.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 1 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_zInt.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 4 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_qzInt.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 5 & 


#################
#Kaon RATIOS
#################

#################
#Rho Ratios
#################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_pip.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 0 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_no_weight.root 1 7 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 0 0 matchCut2D_10.root & 

#################
#Binned RATIOS
#################

#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_pT_qInt.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 5 corrections_10.2_pT_fit.root 1.25 1 0 matchCut2D_10.root 1 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_pT_int.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 1 corrections_10.2_pT_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_pT_no_match.root 0 7 1 acceptanceMap_allE_final.root pT 5 0 1.5 1 corrections_10.2_pT_no_match_fit.root 1.25 1 0 matchCut2D_10.root & 


#################
#RGA RATIOS
#################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_rga_pT_qInt.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 5 corrections_10.2_rga_pT_fit.root 1.25 1 1 matchCut2D_10.root 1 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_rga_pT_int.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 1 corrections_10.2_rga_pT_fit.root 1.25 1 1 matchCut2D_10.root & 

##################
#dOu RATIOS
##################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_rga_qzInt.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 5 corrections_10.2_rga_fit.root 1.25 1 1 matchCut2D_10.root 5 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_qzInt.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 5 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 5 & 


#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_pT.root 1 15 1 acceptanceMap_allE_final.root pT 5 0 1.5 5 corrections_10.2_pT_fit.root 1.25 1 0 matchCut2D_10.root & 

#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_W2_Q2.root 1 15 1 acceptanceMap_allE_final.root W2 7 6 13 7 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 9 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_W2_Mx.root 1 15 1 acceptanceMap_allE_final.root Mx 8 1.7 3.7 8 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 8 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_no_k.root 1 11 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 8 & 

##################
#Corrections RATIOS
###################
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_discMC.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN.root 1.25 1 0 matchCut2D_10.root & 

##################
# RGA RATIOS
##################
#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/final_note/ratios_2d_rho_rga.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_rga_fit.root 1.25 1 1 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/final_note/ratios_2d_rho_rga_qInt.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_rga_fit.root 1.25 1 1 matchCut2D_10.root 1 & 

#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/final_note/ratios_2d_rho_qInt.root 1 15 1 acceptanceMap_allE_final.root pT 1 0 1.6 1 corrections_10.2_AN_fit.root 1.25 1 0 matchCut2D_10.root 1 & 

#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/analysis_note/ratio_rga_rho_pT_qzInt.root 1 15 1 acceptanceMap_allE_final.root pT 4 0 1.2 4 corrections_10.2_rga_pT_fit.root 1.25 1 1 matchCut2D_10.root 5 & 
#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/analysis_note/ratio_10.2_rho_pT_qzInt.root 1 15 1 acceptanceMap_allE_final.root pT 4 0 1.2 4 corrections_10.2_pT_fit.root 1.25 1 0 matchCut2D_10.root 5 & 
#../build/analysis/makeRatio3D_fullUnc 0 ../histograms/analysis_note/ratio_2d_rho_pT.root 1 15 1 acceptanceMap_allE_final.root pT 4 0 1.2 4 corrections_10.2_pT_fit.root 1.25 1 0 matchCut2D_10.root & 
#../build/analysis/makeRatio3D_fullUnc 10.2 ../histograms/analysis_note/ratio_rga_rho_pT.root 1 15 1 acceptanceMap_allE_final.root pT 4 0 1.2 4 corrections_10.2_rga_pT_fit.root 1.25 1 1 matchCut2D_10.root & 