#!/bash

#./skimmers/detector_skimmer 0 10.2 0 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim &
./../build/skimmers/detector_skimmer 3 10.4 5 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim_outbending 0 
#./skimmers/detector_skimmer 0 10.6 0 0 0 /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim &

wait

#./skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root 0 0 10.2 & 
#./../build/skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim_outbending /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim_outbending.root 0 5 10.4 0 
#./skimmers/final_skimmer /volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim /volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root 0 0 10.6 &
