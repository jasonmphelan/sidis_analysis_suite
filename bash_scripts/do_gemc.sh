#!/bash

#./skimmers/detector_skimmer 0 10.2 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim &
#./skimmers/detector_skimmer 0 10.4 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.4/detector_skims/detector_skim &
#./skimmers/detector_skimmer 0 10.6 1 1 0 /volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.6/detector_skims/detector_skim &

wait

./../build/skimmers/final_skimmer ../trees/final_skims/GEMC/detector_skim_10.2.root ../trees/final_skims/GEMC/final_skim_10.2.root 1 1 10.2 &
./../build/skimmers/final_skimmer ../trees/final_skims/GEMC/detector_skim_10.4.root ../trees/final_skims/GEMC/final_skim_10.4.root 1 1 10.4 &
./../build/skimmers/final_skimmer ../trees/final_skims/GEMC/detector_skim_10.6.root ../trees/final_skims/GEMC/final_skim_10.6.root 1 1 10.6 &

