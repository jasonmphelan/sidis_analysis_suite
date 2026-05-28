#!/bash

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.2_mc.root 0 1 corrections_10.2_AN_fit.root null 0 0 1 ../trees/final_skims/10.2/final_skim.root
./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.2_pT_no_match.root 0 1 corrections_10.2_pT_no_match_fit.root pT 4 0 1.2 ../trees/final_skims/10.2/final_skim.root 
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.2_pT.root 1 1 corrections_10.2_pT_fit.root pT 4 0 1.2 ../trees/final_skims/10.2/final_skim.root 

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.4_mc.root 0 1 corrections_10.4_AN_fit.root null 0 0 1 ../trees/final_skims/10.4/final_skim.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.6_mc.root 0 1 corrections_10.6_AN_fit.root null 0 0 1 ../trees/final_skims/10.6/final_skim.root &


#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_out.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/outbending/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_out_etrig.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/outbending/final_skim_etrig.root &

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.2.root 0 0 ../trees/final_skims/GEMC/final_skim_10.2.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.4.root 0 0 ../trees/final_skims/GEMC/final_skim_10.4.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.6.root 0 0 ../trees/final_skims/GEMC/final_skim_10.6.root &

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.2.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10.2.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.4.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10.4.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.6.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10.6.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10637.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10637.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10820.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10820.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10892.root 0 0 1 ../trees/final_skims/GEMC/final_skim_10892.root &


#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10.2.root 0 0 ../trees/final_skims/GEMC/gen_skim_10.2.root &
#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10.4.root 0 0 ../trees/final_skims/GEMC/gen_skim_10.4.root &
#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10.6.root 0 0 ../trees/final_skims/GEMC/gen_skim_10.6.root &
#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10637.root 0 0 ../trees/final_skims/GEMC/gen_skim_10637.root &
#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10820.root 0 0 ../trees/final_skims/GEMC/gen_skim_10820.root &
#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10892.root 0 0 ../trees/final_skims/GEMC/gen_skim_10892.root &


wait

#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_data ../histograms/analysis_note/kinematic_plots_data

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.2.root 0 0 ../trees/final_skims/GEMC/tight_skim_10.2.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.4.root 0 0 ../trees/final_skims/GEMC/tight_skim_10.4.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.6.root 0 0 ../trees/final_skims/GEMC/tight_skim_10.6.root &

#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10.2.root 0 0 ../trees/final_skims/GEMC/gen_skim_comb_10.2.root &
wait
#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_gemc ../histograms/analysis_note/kinematic_plots_gemc
#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_data ../histograms/analysis_note/kinematic_plots_data

wait 

#python ../plotting/plotCharge.py   -o ../plotting/analysis_note/charge_plots '../histograms/analysis_note/kinematic_plots_data_allE.root'
#python ../plotting/plotCharge.py ../histograms/analysis_note/kinematic_plots_data_10.6.root ../plotting/analysis_note/charge_plots/10.6

echo "finished charge"
#python ../plotting/plotEnergy.py ../histograms/analysis_note/kinematic_plots_data ../plotting/analysis_note/energy_kinematic_dep
#python ../plotting/plotOutbending.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data /volatile/clas12/users/jphelan/SIDIS/analysis_note/outbending_kinematics

#python ../plotting/plotDataSim.py ../histograms/analysis_note/kinematic_plots_data_allE.root ../histograms/analysis_note/kinematic_plots_gemc_allE.root ../plotting/analysis_note/data_sim &


#for energy in 10.2 10.4 10.6;
#do
#	python ../plotting/plotDataSim.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_$energy.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_$energy.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/data_sim_plots_$energy &
#done
