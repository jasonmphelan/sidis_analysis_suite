#!/bash

#python ../plotting/plotFragBinned.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_data_binned_q2 ../histograms/analysis_note/ratios_2d_mc.root

#wait

#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_mc_effect ../histograms/analysis_note/ratios_2d_no_corr.root Data  ../histograms/analysis_note/ratios_2d_bin.root '+ Bin Migration' ../histograms/analysis_note/ratios_2d_mc.root '+ Acceptance' &
#python ../plotting/plotRatios.py  /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_k_cut ../histograms/analysis_note/ratios_2d_k.root '$2.5\sigma$' ../histograms/analysis_note/ratios_2d_k_2sig.root '$2\sigma$'  ../histograms/analysis_note/ratios_2d_k_3sig.root '$3\sigma$'
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_k_effect  ../histograms/analysis_note/ratios_2d_mc.root 'MC Corrections' ../histograms/analysis_note/ratios_2d_pi2k.root '+ EB Pion Corrections'  ../histograms/analysis_note/ratios_2d_k.root '+ EB Kaon Corrections'

python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_rho_effect  ../histograms/analysis_note/ratios_2d_k_map.root '+ Kaon Corrections' ../histograms/analysis_note/ratios_2d_rho_no_match_map.root '+ Rho Corrections' &
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_2d_map  ../histograms/analysis_note/ratios_2d_mc_map.root 'acc map + mc' ../histograms/analysis_note/ratios_2d_mc_none.root 'no map + mc' 
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_3d_map  ../histograms/analysis_note/ratios_3d_map.root 'acc map +3D matching(New)' ../histograms/analysis_note/ratios_3d_mc.root 'no map + 3D + mc' ../histograms/analysis_note/ratios_3d_no_corr.root 'no map + 3D (old)'
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_rho_effect  ../histograms/analysis_note/ratios_2d_rho_no_match_map.root '+ No acc matching on rho, no max weight cut' ../histograms/analysis_note/ratios_2d_rho_match_map.root '+ Rho Corrections' ../histograms/analysis_note/ratios_2d_rho_no_max_map.root '+ Rho Corrections without max weight' &


#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_3d_effect ../histograms/analysis_note/ratios_2d_mc.root '2D Matching' ../histograms/analysis_note/ratios_3d_mc.root '3D Matching' &

#wait

#python ../plotting/plotFragBinned3d.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_sector_e ../histograms/analysis_note/ratios_2d_rho_sector_e.root sector_e & 
#python ../plotting/plotFragBinned3d.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_sector_pi ../histograms/analysis_note/ratios_2d_rho_sector_pi.root sector_pi &
#python ../plotting/plotFragBinned3d.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_pT ../histograms/analysis_note/ratios_2d_rho_pT.root pT &

#wait

#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_energy_effect ../histograms/analysis_note/ratios_2d_rho_10.2.root '$E = 10.2$ GeV' ../histograms/analysis_note/ratios_2d_rho_10.4.root '$E = 10.4$ GeV' ../histograms/analysis_note/ratios_2d_rho_10.6.root '$E = 10.6$ GeV'  



