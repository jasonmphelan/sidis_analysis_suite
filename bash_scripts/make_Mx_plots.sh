#!/bash
energy=10.2

#./cutAnalysis/makeMxPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/mx_plots_10.2.root 2 0 10.2
#./cutAnalysis/makeMxPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/mx_plots_10.4.root 2 0 10.4 &
#./cutAnalysis/makeMxPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/mx_plots_10.6.root 2 0 10.6 &

wait

#./analysis/mergeEnergy /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/mx_plots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/mx_plots
python ../plotting/plotMxDep.py ../histograms/analysis_note/mx_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots
