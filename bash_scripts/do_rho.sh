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

for E in 10.6; do
for ((i=27;i<28; i++)); do
#for i in 45 47 48; do
./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_$E.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_$E $E 50 0 acceptanceMap_allE.root $i
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4.root 10.4 10 0 acceptanceMap_allE.root &
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6.root 10.6 10 0 acceptanceMap_allE.root &

done
wait
done
wait

#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2_loose.root 10.2 10 0 acceptanceMap_allE_loose.root &
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4_loose.root 10.4 10 0 acceptanceMap_allE_loose.root &
#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6_loose.root 10.6 10 0 acceptanceMap_allE_loose.root &


wait

#./rhoAnalysis/norms _acc /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms_acc &
#./rhoAnalysis/norms "" /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms &
#./rhoAnalysis/norms _loose /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms_loose &

wait

#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back.root &
#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back_acc.root _acc &
#./rhoAnalysis/studyBack /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_back_loose.root _loose &



wait

#python ../plotting/plotRhoBac.py ../histograms/analysis_note/rho_back.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_background &
#python ../plotting/plotRhoFits.py ../histograms/analysis_note/rho_fits.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_scaling &
#llpython ../plotting/plotRhoSub.py ../histograms/analysis_note/rho_back.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_subtraction &

#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.2.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2_acc.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.4.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4_acc.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_10.6.root 1 0 /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6_acc.root &

#wait

#./analysis/mergeEnergy /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho
#python ../plotting/plotCharge.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_rho_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/rho_kinematics
