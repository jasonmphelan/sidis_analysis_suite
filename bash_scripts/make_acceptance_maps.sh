#!/bash

./cutAnalysis/makeThetaPhi /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.2.root 0 0 10.2 &
./cutAnalysis/makeThetaPhi /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.4.root 0 0 10.4 &
./cutAnalysis/makeThetaPhi /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.6.root 0 0 10.6 &

wait

#./cutAnalysis/makeAcceptanceMap /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.2.root /work/clas12/users/jphelan/sidis_analysis_suite/data/acceptance_map/acceptanceMap_10.2.root
#./cutAnalysis/makeAcceptanceMap /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.4.root /work/clas12/users/jphelan/sidis_analysis_suite/data/acceptance_map/acceptanceMap_10.4.root
#./cutAnalysis/makeAcceptanceMap /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/theta_phi_10.6.root /work/clas12/users/jphelan/sidis_analysis_suite/data/acceptance_map/acceptanceMap_10.6.root
