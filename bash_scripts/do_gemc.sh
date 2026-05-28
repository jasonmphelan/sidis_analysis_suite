#!/bash

#./../build/skimmers/detector_skimmer 0 10.2 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim 1 &
#./../build/skimmers/detector_skimmer 10 10.4 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/detector_skims/detector_skim_test 0 &
#./../build/skimmers/detector_skimmer 0 10.6 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/detector_skims/detector_skim_test 0

wait

./../build/skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim.root /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/final_skims/inal_skim.root 1 1 10.2 1
#./../build/skimmers/final_skimmer ../trees/final_skims/GEMC/detector_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10.4.root 1 1 10.4 
#./../build/skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/detector_skims/detector_skim_test.root  /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/final_skims/final_skim.root 1 1 10.6 0 

