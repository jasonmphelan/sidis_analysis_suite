import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys
plt.rcParams["font.family"] = "sans-serif"
colorList = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']
colorList = ['r','b','o']

def setTitle2D( ax, key):
	if "Beta" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]")
		ax.set_ylabel(r"$\beta$")
	if "WV" in key:
		ax.set_xlabel(r"V$_{PCAL}$ [cm]")
		ax.set_ylabel(r"W$_{PCAL}$ [cm]")
	if "Edep" in key:
		ax.set_xlabel(r"$E_{PCAL}$ [GeV]")
		ax.set_ylabel(r"E_{CIN} + E_{COUT} [GeV]")
	if "SF_corr" in key:
		ax.set_xlabel(r"$\delta E_{dep}(PCAL)/p_{e}$")
		ax.set_ylabel(r"$(E_{PCAL} + E_{CIN} + E_{COUT})/p_e$ [GeV]")
	if "SF_sec" in key:
		ax.set_xlabel(r"$p_{e}$ [GeV]")
		ax.set_ylabel(r"$\delta E_{CIN}/p_e$")
	if "Chi2" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]")
		ax.set_ylabel(r"$\chi^2_{PID}$")
	if "Chi2" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]")
		ax.set_ylabel(r"$\chi^2_{PID}$")
	if "Q2_Z" in key:
		ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
		ax.set_ylabel(r"$z$")
	if "Q2_W" in key:
		ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
		ax.set_ylabel(r"$W$ [GeV]")
	if "Q2_omega" in key:
		ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
		ax.set_ylabel(r"$\omega$ [GeV]")
def setTitle1D( ax, key):
	ax.axhline(y=1, linestyle='--',c='k', linewidth=2)
	if "Beta" in key:
		ax.set_xlabel(r"$\beta$", fontsize=16)
	if "Vz_e" in key:
		ax.set_xlabel(r"$v_z^e$ [cm]", fontsize=16)
	if "Vz_pi" in key:
		ax.set_xlabel(r"$v_z^{\pi} - v_z^e$ [cm]", fontsize=16)
	if "Z_pi" in key:	
		ax.set_xlabel(r"$z$", fontsize=16)
	if "Y_pi" in key:
		ax.set_xlabel(r"$y$", fontsize=16)
	if "W_pi" in key:
		ax.set_xlabel(r"$W$ [GeV]", fontsize=16)
	if "Theta_e" in key:
		ax.set_xlabel(r"$\theta_e$ [deg.]", fontsize=16)
	if "Theta_pi" in key:
		ax.set_xlabel(r"$\theta_{\pi} [deg.]$", fontsize=16)
	if "Q2" in key:
		ax.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=16)
	if "P_pi" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=16)
	if "P_e" in key:
		ax.set_xlabel(r"$p_{e}$ [GeV]", fontsize=16)
	if "Mx" in key:
		ax.set_xlabel(r"$M_X$ [GeV]", fontsize=16)
	if "Omega" in key:
		ax.set_xlabel(r"$\omega$ [GeV]", fontsize=16)

#inFile_names = ['../histograms/analysis_note/kinematics_data_10.4_p_cut.root', '../histograms/analysis_note/kinematics_kaons_10.2_p_cut.root', '../histograms/analysis_note/kinematics_data_10.4_beta_cut.root', '../histograms/analysis_note/kinematics_kaons_10.2_beta_cut.root']

inFile_names = sys.argv[1]
out_names = sys.argv[2]

comp_files = {0,1}
#inFile_names = ['../histograms/analysis_note/kinematics_data_10.4_p_cut.root', '../histograms/analysis_note/kinematics_data_10.4_beta_cut.root']
labels = [r'$Y_{\pi}^{EB}$', r'$Y_{k}^{EB}$',r'$Y_{\pi}^{EB}$, Hit in RICH', r'$Y_{k}^{EB}$, Hit in RICH']

inFile = uproot.open(inFile_names) 

for key in inFile.keys():
	if "pim" in key:
		continue
	if "bac" not in key:
		continue
	if isinstance( inFile[key], uproot.models.TH.Model_TH1F_v3): 
		
		pip = inFile[key]
		pip_sig = inFile[ key.replace("_bac", "")]
		hName = key.replace("_bac", "")
		pim = inFile[key.replace("pip", "pim")]
		pim_sig = inFile[ hName.replace("pip", "pim")]
		
		values_pip = pip.values()		
		errors_pip = pip.errors()
		values_pip_sig = pip_sig.values()		
		errors_pip_sig = pip_sig.errors()

		values_pim = pim.values()
		errors_pim = pim.errors()
		values_pim_sig = pim_sig.values()		
		errors_pim_sig = pim_sig.errors()

		if np.sum( values_pim ) < 10 or np.sum(values_pip) < 10:
			continue

		#values_pip[values_pip == 0] = np.nan
		#values_pim[values_pim == 0] = np.nan
		#values_pip_sig[values_pip_sig == 0] = np.nan
		#values_pim_sig[values_pim_sig == 0] = np.nan
		
		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2

		#values_pip = (values_pip_sig[np.where(binCenters==.45)[0]]/values_pip[np.where(binCenters==.45)[0]])*values_pip	
		#values_pim = (values_pim_sig[np.where(binCenters==.45)[0]]/values_pim[np.where(binCenters==.45)[0]])*values_pim	

		indices = np.where( binCenters < 0.525)

		scale_pip = np.sum(values_pip_sig[indices[0]])/np.sum(values_pip[indices[0]])
		scale_pim = np.sum(values_pim_sig[indices[0]])/np.sum(values_pim[indices[0]])
		

		fig, ax = plt.subplots(2, figsize=(12,6), layout='constrained', sharex=True)

		pip_sub = values_pip_sig - values_pip*scale_pip
		pip_sub[pip_sub <= 0] = np.nan

		pim_sub = values_pim_sig - values_pim*scale_pim
		pim_sub[pim_sub <= 0] = np.nan

		ax[0].plot(binCenters, pip_sub, color='b',drawstyle='steps-mid', label = r'$\pi^+$')
		ax[0].plot(binCenters, pim_sub, color='r',drawstyle='steps-mid', label = r'$\pi^-$')

		ratio = (pip_sub)/ (pim_sub)
	


		#ax[0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
		ax[0].set_ylabel('Counts [a.u.]', fontsize=16)
		ax[1].set_ylabel('Counts [a.u.]', fontsize=16)
		ax[1].set_xlabel(r'$M^{\pi+\pi-}$ [GeV]', fontsize=16)
		ax[0].set_title(r'$(e, e^{\prime}\pi^+)$', fontsize=16)	
		ax[1].set_title(r'$\pi^-$', fontsize=16)	

		ax[1].plot(binCenters, ratio,  color='k', linestyle=' ',marker='o')
		ax[1].axhline(y=1, linestyle='--',c='k', linewidth=2)
		ax[1].set_ylim(0.5, 2)
		#ax[1].set_ylim(0.5, 1.5)
		#ax[1].plot(binCenters, values_pim_sub, color='k',drawstyle='steps-mid')#, label = r'$\pi^-$')

		#ax[1].plot(binCenters, values_pim,  color='r', drawstyle='steps-mid')

			#ax[1, 1].errorbar(binCenters, ratio_pim, ratio_err_pim, capsize=2, fmt='k o', ecolor='k')
			#setTitle1D( ax[1,1], key )
			
			

		plt.margins(y=0)
		plt.margins(x=0)
			
		for j in range(2):
			ax[j].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			ax[j].tick_params(which='both', width=2, labelsize=12)
			ax[j].tick_params(which='major', length=7, labelsize=12)
			ax[j].tick_params(which='minor', length=4, labelsize=12)
			
			#ax[0].set_ylabel("Counts [a.u.]", fontsize=18)
			#ax[1].set_ylabel(r"$Y(e, e'\pi^+)/Y(e, e'\pi^-)$", fontsize=12)
			#setTitle1D( ax[1], key )

	pdfName = hName.replace("_pip", "")
	fig.legend(fontsize=16)
	fig.savefig(f"{out_names}/{pdfName}.pdf")
	plt.close(fig)