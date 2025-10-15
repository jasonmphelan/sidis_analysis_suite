#!/bash
energy=10.2

#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.2 10.2 0 &
#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.4 10.4 0 &
#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.6 10.6 0 &
./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_out_10.4 10.4 5 &
wait
#./analysis/mergeEnergy /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots
#python ../plotting/plotCuts.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_allE
#python ../plotting/plotSF.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_allE
#python ../plotting/plotFiducials.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_allE

python ../plotting/plotCuts.py ../histograms/analysis_note/selection_out_10.4.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_out
python ../plotting/plotSF.py ../histograms/analysis_note/selection_out_10.4.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_out
python ../plotting/plotFiducials.py ../histograms/analysis_note/selection_out_10.4.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_out
