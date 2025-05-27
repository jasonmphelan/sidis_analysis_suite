#!/bash
energy=10.2

./cutAnalysis/makeKinematicPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_cuts_10.2 0 0 10.2 &
./cutAnalysis/makeKinematicPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_cuts_10.4 0 0 10.4 &
./cutAnalysis/makeKinematicPlots /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_cuts_10.6 0 0 10.6 &

wait

./analysis/mergeEnergy ../histograms/analysis_note/kinematic_cuts ../histograms/analysis_note/kinematic_cuts

python ../plotting/plotCuts.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/kinematic_cuts_allE.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/kinematic_cuts
