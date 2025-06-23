#!/bash

#./skimmers/detector_skimmer 0 10.2 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim &
./skimmers/detector_skimmer 0 10.4 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/detector_skims/detector_skim 
#./skimmers/detector_skimmer 0 10.6 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/detector_skims/detector_skim &

wait

#./skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim.root /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/final_skims/final_skim.root 1 1 10.2 &
./skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/detector_skims/detector_skim.root /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/final_skims/final_skim.root 1 1 10.4 
#./skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/detector_skims/detector_skim.root /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root 1 1 10.6 &

