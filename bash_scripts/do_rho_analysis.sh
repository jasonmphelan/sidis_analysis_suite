#!/bash

#./skimmers/detector_skimmer 0 10.4 4 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/outbending/run_skim
#./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root 0 0 10.2

#./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/final_skim_10.2.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2.root 10.2 1
./rhoAnalysis/norms x /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms.cpp
