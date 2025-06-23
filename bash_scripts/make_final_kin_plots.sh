#!/bash

#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_10.2.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_10.4.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_10.6.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_out.root 0 0 /volatile/clas12/users/jphelan/SIDIS/data/final_skims/outbending/final_skim.root &

wait

#./analysis/mergeEnergy /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data

#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_10.2.root 0 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/final_skims/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_10.4.root 0 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/final_skims/final_skim.root &
#./analysis/makeKinematics /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_10.6.root 0 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root &

wait

#python ../plotting/plotCharge.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/charge_plots

python ../plotting/plotEnergy.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data /volatile/clas12/users/jphelan/SIDIS/analysis_note/energy_kinematic_dep

#for energy in 10.2 10.4 10.6;
#do
#	python ../plotting/plotDataSim.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_data_$energy.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_plots_gemc_$energy.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/data_sim_plots_$energy
#done
