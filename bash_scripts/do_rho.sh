#!/bash

#./skimmers/detector_skimmer 0 10.4 4 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/outbending/run_skim
#./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root 0 0 10.2 &
#./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root 0 0 10.4 &
#./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root 0 0 10.6 &

wait

#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2_acc.root 10.2 10 1 acceptanceMap_allE.root & 
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4_acc.root 10.4 10 1 acceptanceMap_allE.root &
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6_acc.root 10.6 10 1 acceptanceMap_allE.root &

#wait

for E in 10.2; do
for ((i=0;i<1; i++)); do
#for i in 45 47 48; do
./rhoAnalysis/rotateRhoSym /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_$E.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_sym_$E $E 20 1 acceptanceMap_allE_final.root $i 

done
#wait
done
#wait


#for bound in 1.3 1.4 1.5 1.6 2.5; do
#(
#    ./rhoAnalysis/norms "_acc" /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms_Mx_$bound 1.15 $bound 0.525
#    ./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_match_map_Mx_${bound}.root 2 5 1 acceptanceMap_allE_final.root rho_norms_Mx_${bound}_acc.root 
#)&
#done
#wait

#./rhoAnalysis/norms "_acc" /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms_no_weight 1.15 1.6 0.525
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_match_map_no_weight.root 2 5 1 acceptanceMap_allE_final.root rho_norms_no_weight_acc.root 

#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_rho_norm_cut_effect  ../histograms/analysis_note/ratios_2d_rho_match_map_0.475.root  'Norm Bin = 0.475' ../histograms/analysis_note/ratios_2d_rho_match_map_0.500.root  'Norm Bin = 0.500' ../histograms/analysis_note/ratios_2d_rho_match_map_0.525.root  'Norm Bin = 0.525' ../histograms/analysis_note/ratios_2d_rho_match_map_0.550.root  'Norm Bin = 0.550'
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_rho_norm_bounds_effect  ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_1.3.root  'Mx max = 1.3' ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_1.4.root  'Mx max = 1.4' ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_1.5.root  'Mx max = 1.5' ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_1.6.root  'Mx max = 1.6' ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_2.5.root  'Mx max = 2.5'
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_rho_weight_effect  ../histograms/analysis_note/ratios_2d_rho_match_map_Mx_1.6.root  'With Weights' ../histograms/analysis_note/ratios_2d_rho_match_map_no_weight.root  'Without Weights'
#./rhoAnalysis/norms "" /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms &
#./rhoAnalysis/norms _loose /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms_loose &

wait

#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back.root 
#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back_acc.root _acc &
#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back_loose.root _loose &



wait

#python ../plotting/plotRhoBac.py ../histograms/analysis_note/rho_back_acc.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_background &
#python ../plotting/plotRhoFits.py ../data/correctionFiles/rho_norms_M_rho.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_scaling &
#python ../plotting/plotRhoSub.py ../histograms/analysis_note/rho_back_acc.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_subtraction 

#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.2.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2_acc.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.4.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4_acc.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.6.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6_acc.root &

#wait

#./analysis/mergeEnergy /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho
#python ../plotting/plotCharge.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_kinematics
