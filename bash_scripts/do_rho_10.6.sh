#!/bash

./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root 0 0 10.6

./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/final_skim_10.6.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.26.root 10.6 1
