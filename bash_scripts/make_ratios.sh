#!/bash

#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_raw.root 0 0 &
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_no_corr.root 2 0 &
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_3d_no_corr.root 3 0 &

#wait

#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_bin.root 2 1 &
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_acc.root 2 2 &

#wait

#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_mc.root 2 3 
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_3d_mc.root 3 3 &
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_k.root 2 4 
#./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_k_2sig.root 2 4 

wait

./analysis/makeRatio 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho.root 2 5 &
#./analysis/makeRatio 10.2 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_10.2.root 2 5 & 
#./analysis/makeRatio 10.4 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_10.4.root 2 5 &
#./analysis/makeRatio 10.6 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_10.6.root 2 5 &

wait

#./analysis/makeRatioBinned 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_pT.root 2 5 pT 
#./analysis/makeRatioBinned 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_sector_e.root 2 5 sector_e
#./analysis/makeRatioBinned 0 /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/ratios_2d_rho_sector_pi.root 2 5 sector_pi
