#!/bash

#python ../plotting/plotFragBinned.py ../plotting/analysis_note/ratio_data_binned_q2 ../histograms/analysis_note/ratios_2d_rho_pim.root &
#python ../plotting/plotRatios_xB.py ../plotting/analysis_note/ratio_data_binned_xB ../histograms/analysis_note/ratios_2d_rho_pim.root &

#python ../plotting/plotFragBinned.py ../plotting/analysis_note/ratio_data_binned_raw ../histograms/analysis_note/ratios_2d_map_no_corr.root
#python ../plotting/plotFragBinned.py ../plotting/analysis_note/ratio_data_binned_mc ../histograms/analysis_note/ratios_2d_acc.root
#python ../plotting/plotFragBinned.py ../plotting/analysis_note/ratio_data_binned_k ../histograms/analysis_note/ratios_2d_k.root


#wait

#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_mc_effect ../histograms/analysis_note/ratios_2d_no_corr.root Data  ../histograms/analysis_note/ratios_2d_bin.root '+ Bin Migration' ../histograms/analysis_note/ratios_2d_mc.root '+ Acceptance' &
#python ../plotting/plotRatios.py  ../plotting/analysis_note/ratio_k_sigma ../histograms/analysis_note/ratios_2d_k.root '$2.5\sigma$' ../histograms/analysis_note/ratios_2d_k_2sig.root '$2\sigma$' ../histograms/analysis_note/ratios_2d_k_3sig.root '$3\sigma$' &
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_k_effect  ../histograms/analysis_note/ratios_2d_mc.root 'MC Corrections' ../histograms/analysis_note/ratios_2d_pi2k.root '+ EB Pion Corrections'  ../histograms/analysis_note/ratios_2d_k.root '+ EB Kaon Corrections'

#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_effect  ../histograms/analysis_note/ratios_2d_k_match_map.root '+ Kaon Corrections' ../histograms/analysis_note/ratios_2d_rho_match_map_no_fix.root '+ Rho Corrections' ../histograms/analysis_note/ratios_2d_rho_match_map_fix.root '+ Rho Corrections Fixed' 
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_2d_map  ../histograms/analysis_note/ratios_2d_mc_map.root 'acc map + mc' ../histograms/analysis_note/ratios_2d_mc_none.root 'no map + mc' 
#python ../plotting/plotRatios.py /volatile/clas12/users/jphelan/SIDIS/analysis_note/ratio_3d_map  ../histograms/analysis_note/ratios_3d_map.root 'acc map +3D matching(New)' ../histograms/analysis_note/ratios_3d_mc.root 'no map + 3D + mc' ../histograms/analysis_note/ratios_3d_no_corr.root 'no map + 3D (old)'
#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_effect  ../histograms/analysis_note/ratios_2d_k.root '+ k' ../histograms/analysis_note/ratios_2d_rho_no_sym.root '+ No sym on rho' ../histograms/analysis_note/ratios_2d_rho_no_sym_test.root '+ Rho Corrections' &

#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_effect ../histograms/analysis_note/ratios_2d_rho_no_sub.root '+ VM (no sub)' ../histograms/analysis_note/ratios_2d_rho_pim.root '+ VM (pim)' ../histograms/analysis_note/ratios_2d_rho_pip.root '+ VM (pip)'
#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_effect  ../histograms/analysis_note/ratios_2d_k.root 'Data + MC Corrections + K Corrections' ../histograms/analysis_note/ratios_2d_rho_pim.root '+ VM Subtraction'




#wait

#python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_sector_e ../histograms/analysis_note/ratios_2d_rho_sector_e.root sector_e & 
#python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_sector_pi ../histograms/analysis_note/ratios_2d_rho_sector_pi.root sector_pi & 
#python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_pT ../histograms/analysis_note/ratios_2d_rho_pt.root pT & 

python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_phi_q ../histograms/analysis_note/ratios_2d_phi_q.root phi_q & 
#python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_eta ../histograms/analysis_note/ratios_2d_rho_eta.root eta & 


python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_sector_pi ../histograms/analysis_note/ratios_2d_sector_pi.root sector_pi &
python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_sector_e ../histograms/analysis_note/ratios_2d_sector_e.root sector_e &

python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_pT ../histograms/analysis_note/ratios_2d_pT.root pT &

#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_sector_pi ../histograms/analysis_note/ratios_2d_acc.root '$E = 10.2$ GeV' ../histograms/analysis_note/ratios_no_match_mc.root 'With sector corr' 

#wait
#
#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_energy_effect ../histograms/analysis_note/ratios_2d_rho_10.2.root '$E = 10.2$ GeV' ../histograms/analysis_note/ratios_2d_rho_10.4.root '$E = 10.4$ GeV' ../histograms/analysis_note/ratios_2d_rho_10.6.root '$E = 10.6$ GeV'  &


#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_mx_effect  ../histograms/analysis_note/ratios_2d_acc.root '$M_X > 1.7$ [GeV]' ../histograms/analysis_note/ratios_mid_mx.root '$M_X > 1.6$ [GeV]'  ../histograms/analysis_note/ratios_loose_mx.root '$M_X > 1.5$ [GeV]' &
python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_cut_effect  ../histograms/analysis_note/ratios_2d_rho.root '$M_X^{\pi+, \pi-} < 1.25$ [GeV]' ../histograms/analysis_note/ratios_2d_mx_1.1.root '$M_X^{\pi+, \pi-} < 1.1$ [GeV]' ../histograms/analysis_note/ratios_2d_mx_1.2.root '$M_X^{\pi+, \pi-} < 1.2$ [GeV]' ../histograms/analysis_note/ratios_2d_mx_1.3.root '$M_X^{\pi+, \pi-} < 1.3$ [GeV]' ../histograms/analysis_note/ratios_2d_mx_1.4.root '$M_X^{\pi+, \pi-} < 1.4$ [GeV]' & 

#python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_acc_effect  ../histograms/analysis_note/ratios_2d_rho.root 'Old Code$ [GeV]' ../histograms/analysis_note/ratios_2d_rho_pim.root 'New Code'

python ../plotting/plotRatios.py ../plotting/analysis_note/ratio_rho_acc_effect  ../histograms/analysis_note/ratios_2d_rho_pim.root 'Tagging $\pi^+$' ../histograms/analysis_note/ratios_2d_rho_pip.root 'Tagging $\pi^-$' 

python ../plotting/plotFragBinned3d.py ../plotting/analysis_note/ratio_W2 ../histograms/analysis_note/ratios_2d_rho_W2.root W2 &
