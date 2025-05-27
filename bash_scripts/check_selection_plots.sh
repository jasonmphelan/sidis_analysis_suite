#!/bash
energy=10.2

./cutAnalysis/checkSelection /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/test_selection_plots_$1 $1 1
python ../plotting/plotCuts.py ../histograms/analysis_note/test_selection_plots_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/check_selection
#python ../plotting/plotSF.py ../histograms/analysis_note/selection_plots_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_$1
#python ../plotting/plotFiducials.py ../histograms/analysis_note/selection_plots_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_$1
