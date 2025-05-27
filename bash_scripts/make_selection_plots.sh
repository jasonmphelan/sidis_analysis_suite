#!/bash
energy=10.2

#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.2 10.2 0 &
#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.4 10.4 0 &
#./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_plots_10.6 10.6 0 &

#wait

#./analysis/mergeEnergy ../histograms/analysis_note/selection_plots ../histograms/analysis_note/selection_plots

python ../plotting/plotCuts.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots
python ../plotting/plotSF.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots pip
python ../plotting/plotSF.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots pim
python ../plotting/plotFiducials.py ../histograms/analysis_note/selection_plots_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots
