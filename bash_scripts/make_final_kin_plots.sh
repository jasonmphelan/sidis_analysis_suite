#!/bash

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.2.root 0 0 ../trees/final_skims/10.2/final_skim.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.4.root 0 0 ../trees/final_skims/10.4/final_skim.root &
./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_data_10.6.root 0 0 ../trees/final_skims/10.6/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_out.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/outbending/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_out_etrig.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/outbending/final_skim_etrig.root &



#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_data ../histograms/analysis_note/kinematic_plots_data

#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.2.root 0 0 ../trees/final_skims/GEMC/final_skim_comb_10.2.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.4.root 0 0 ../trees/final_skims/GEMC/final_skim_10.4.root &
#./../build/analysis/makeKinematics ../histograms/analysis_note/kinematic_plots_gemc_10.6.root 0 0 ../trees/final_skims/GEMC/final_skim_10.6.root &

#./../build/analysis/makeGenKinematics ../histograms/analysis_note/kinematic_plots_gen_10.2.root 0 0 ../trees/final_skims/GEMC/gen_skim_comb_10.2.root &
wait
#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_gemc ../histograms/analysis_note/kinematic_plots_gemc
#./../build/analysis/mergeEnergy ../histograms/analysis_note/kinematic_plots_data ../histograms/analysis_note/kinematic_plots_data

#python ../plotting/plotCharge.py ../histograms/analysis_note/kinematic_plots_data_allE.root ../plotting/analysis_note/charge_plots
python ../plotting/plotCharge.py ../histograms/analysis_note/kinematic_plots_data_10.6.root ../plotting/analysis_note/charge_plots/10.6

#python ../plotting/plotEnergy.py ../histograms/analysis_note/kinematic_plots_data ../plotting/analysis_note/energy_kinematic_dep
#python ../plotting/plotOutbending.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data /volatile/clas12/users/jphelan/SIDIS/analysis_note/outbending_kinematics

#python ../plotting/plotDataSim.py ../histograms/analysis_note/kinematic_plots_data_allE.root ../histograms/analysis_note/kinematic_plots_gemc_allE.root ../plotting/analysis_note/data_sim &


#for energy in 10.2 10.4 10.6;
#do
#	python ../plotting/plotDataSim.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_$energy.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_$energy.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/data_sim_plots_$energy &
#done
