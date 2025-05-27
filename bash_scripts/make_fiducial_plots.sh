#!/bash
energy=10.2

./cutAnalysis/fiducials /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/fiducial_plots_$1.root 0 0 $1
python ../plotting/plotEdge.py ../histograms/analysis_note/fiducial_plots_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_plots_$1
