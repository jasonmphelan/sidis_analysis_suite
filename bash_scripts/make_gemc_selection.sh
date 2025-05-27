#!/bash
energy=10.2

./cutAnalysis/makeDetectorPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/selection_gemc_$1 $1 1500
python ../plotting/plotCuts.py ../histograms/analysis_note/selection_gemc_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_gemc_$1
python ../plotting/plotSF.py ../histograms/analysis_note/selection_gemc_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_gemc_$1
python ../plotting/plotFiducials.py ../histograms/analysis_note/selection_gemc_$1.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/selection_gemc_$1
