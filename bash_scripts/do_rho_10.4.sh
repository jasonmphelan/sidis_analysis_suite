#!/bash

#./skimmers/detector_skimmer 0 10.4 4 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/outbending/run_skim
./skimmers/rho_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root 0 0 10.4

./rhoAnalysis/rotateRho /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/final_skim_10.4.root /volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4.root 10.4 1
