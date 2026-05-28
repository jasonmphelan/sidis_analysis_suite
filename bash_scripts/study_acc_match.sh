#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_1.root 0 .01 &

#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_5.root 0 .05 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_10.root 0 .1 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_15.root 0 .15 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_20.root 0 .20 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_25.root 0 .25 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_50.root 0 .50 &

wait

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_AN_test.root 1 0 0 matchCut2D_1.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_AN_test.root 1 0 0 matchCut2D_1.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_AN_test.root 1 0 0 matchCut2D_1.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_5.root 1 0 0 matchCut2D_5.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_5.root 1 0 0 matchCut2D_5.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_5.root 1 0 0 matchCut2D_5.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_10.root 1 0 0 matchCut2D_10.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_10.root 1 0 0 matchCut2D_10.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_10.root 1 0 0 matchCut2D_10.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_15.root 1 0 0 matchCut2D_15.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_15.root 1 0 0 matchCut2D_15.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_15.root 1 0 0 matchCut2D_15.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_20.root 1 0 0 matchCut2D_20.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_20.root 1 0 0 matchCut2D_20.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_20.root 1 0 0 matchCut2D_20.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_25.root 1 0 0 matchCut2D_25.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_25.root 1 0 0 matchCut2D_25.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_25.root 1 0 0 matchCut2D_25.root 0 &

./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.2.root ../trees/final_skims/GEMC/gen_skim_comb_10.2.root ../data/correctionFiles/corrections_10.2_50.root 1 0 0 matchCut2D_50.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.4.root ../trees/final_skims/GEMC/gen_skim_10.4.root ../data/correctionFiles/corrections_10.4_50.root 1 0 0 matchCut2D_50.root 0 &
./../build/corrections/calcCorrections ../trees/final_skims/GEMC/tight_skim_10.6.root ../trees/final_skims/GEMC/gen_skim_10.6.root ../data/correctionFiles/corrections_10.6_50.root 1 0 0 matchCut2D_50.root 0 &

wait

./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_1.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_AN_test.root 1.25 1 0 matchCut2D_1.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_5.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_5.root 1.25 1 0 matchCut2D_5.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_10.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_10.root 1.25 1 0 matchCut2D_10.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_15.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_15.root 1.25 1 0 matchCut2D_15.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_20.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_20.root 1.25 1 0 matchCut2D_20.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_25.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_25.root 1.25 1 0 matchCut2D_25.root & 
./../build/analysis/makeRatio3D 0 ../histograms/analysis_note/ratios_2d_50.root 2 3 1 acceptanceMap_allE_final.root Mx 1 1.7 3.7 1 corrections_10.2_50.root 1.25 1 0 matchCut2D_50.root &

wait


python3 ../plotting/plotRatios.py ../plotting/analysis_note/ratio_acc_match \
    ../histograms/analysis_note/ratios_2d_1.root   "1%" \
    ../histograms/analysis_note/ratios_2d_5.root   "5%" \
    ../histograms/analysis_note/ratios_2d_10.root  "10%" \
    ../histograms/analysis_note/ratios_2d_15.root  "15%" \
    ../histograms/analysis_note/ratios_2d_20.root  "20%" \
    ../histograms/analysis_note/ratios_2d_25.root  "25%" \
    ../histograms/analysis_note/ratios_2d_50.root  "50%"

