import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
plt.rcParams["font.family"] = "sans-serif"
colorList = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']

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

inFile_names = ['../histograms/analysis_note/kinematics_data_10.4_p_cut.root', '../histograms/analysis_note/kinematics_kaons_10.2_p_cut.root', '../histograms/analysis_note/kinematics_data_10.4_beta_cut.root', '../histograms/analysis_note/kinematics_kaons_10.2_beta_cut.root']

comp_files = {1, 3}
#inFile_names = ['../histograms/analysis_note/kinematics_data_10.4_p_cut.root', '../histograms/analysis_note/kinematics_data_10.4_beta_cut.root']
labels = [r'$Y_{\pi}^{EB}$', r'$Y_{k}^{EB}$',r'$Y_{\pi}^{EB}$, Hit in RICH', r'$Y_{k}^{EB}$, Hit in RICH']

inFile = [ uproot.open(i) for i in inFile_names ]

for key in inFile[0].keys():
	if "pim" in key:
		continue
	
	if isinstance( inFile[0][key], uproot.models.TH.Model_TH1F_v3): 
		fig, ax = plt.subplots(2, 2, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		values_pip = np.empty( len(inFile), dtype='object')
		values_pim = np.empty( len(inFile), dtype='object')
		errors_pip = np.empty( len(inFile), dtype='object')
		errors_pim = np.empty( len(inFile), dtype='object')
		
		for i in range(len(inFile)):
			file = inFile[i]	
			
			pip = file[key]
			hName = key
			pim = file[hName.replace("pip", "pim")]

			values_pip[i] = pip.values()		
			errors_pip[i] = pip.errors()
			
			values_pim[i] = pim.values()
			errors_pim[i] = pim.errors()

			values_pip[i][values_pip[i] == 0] = np.nan
			values_pim[i][values_pim[i] == 0] = np.nan
			
			values_pip[i] = values_pip[i]/np.nansum(pip.values())
			errors_pip[i] = errors_pip[i]/np.nansum(pip.values())
			values_pim[i] = values_pim[i]/np.nansum(pim.values())
			errors_pim[i] = errors_pim[i]/np.nansum(pim.values())

			ratio_list = [0, 0, 2, 2]

			ratio_pip = values_pip[i]/values_pip[ratio_list[i]]
			ratio_pim = values_pim[i]/values_pim[ratio_list[i]]
			
			ratio_err_pip = ratio_pip*np.sqrt( (errors_pip[i]/values_pip[i])**2 + (errors_pip[0]/values_pip[0])**2)
			ratio_err_pim = ratio_pim*np.sqrt( (errors_pim[i]/values_pim[i])**2 + (errors_pim[0]/values_pim[0])**2)

			xEdges = pip.axis().edges()
			binCenters = (xEdges[:-1]+xEdges[1:])/2

			ax[0, 0].errorbar(binCenters, values_pip[i], errors_pip[i], capsize=2, fmt=colorList[i], ecolor=colorList[i],drawstyle='steps-mid', label=labels[i])
			ax[0,0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
			ax[0,0].set_ylabel('Counts [a.u.]', fontsize=16)
			if i == 1 or i == 3:
				ax[1, 0].errorbar(binCenters, ratio_pip, ratio_err_pip, capsize=2, fmt=colorList[i], ecolor=colorList[i], marker='o')
				setTitle1D( ax[1,0], key )	
				ax[1, 0].set_ylabel(r'$Y_k/Y\pi$', fontsize=16)		

			ax[0, 1].errorbar(binCenters, values_pim[i], errors_pim[i], capsize=2, fmt=colorList[i], ecolor=colorList[i],drawstyle='steps-mid')
			ax[0,1].set_title(rf"$d(e, e'\pi^-)$", fontsize=16)
			if i == 1 or i == 3:
				ax[1, 1].errorbar(binCenters, ratio_pim, ratio_err_pim, capsize=2, fmt=colorList[i], ecolor=colorList[i], marker='o')
				setTitle1D( ax[1,1], key )	
				ax[1, 1].set_ylabel(r'$Y_k/Y_{\pi}$', fontsize=16)		
			#setTitle1D( fig, key )
			#ax[1, 1].errorbar(binCenters, ratio_pim, ratio_err_pim, capsize=2, fmt='k o', ecolor='k')
			#setTitle1D( ax[1,1], key )
			
			

			plt.margins(y=0)
			plt.margins(x=0)
			
			for j in range(2):
				for i in range(2):
					ax[j,i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
					ax[j,i].tick_params(which='both', width=2, labelsize=12)
					ax[j,i].tick_params(which='major', length=7, labelsize=12)
					ax[j,i].tick_params(which='minor', length=4, labelsize=12)
			
			#ax[0].set_ylabel("Counts [a.u.]", fontsize=18)
			#ax[1].set_ylabel(r"$Y(e, e'\pi^+)/Y(e, e'\pi^-)$", fontsize=12)
			#setTitle1D( ax[1], key )

		pdfName = hName.replace("_pip", "")
		fig.legend(fontsize=16)
		fig.savefig(f"comp_plots/{pdfName}.pdf")
		plt.close(fig)

	if isinstance( inFile[0][key], uproot.models.TH.Model_TH2F_v4): 
		fig, axs = plt.subplots(2, 2, figsize=(12,6), height_ratios=[1,1], sharex='col', layout='constrained')
		axs[0,0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
		axs[0,1].set_title(rf"$d(e, e'\pi^-)$", fontsize=16)
			
		
		hist_pip = np.empty( len(inFile), dtype='object')
		hist_pim = np.empty( len(inFile), dtype='object')
		values_pip = np.empty( len(inFile), dtype='object')
		values_pim = np.empty( len(inFile), dtype='object')
		errors_pip = np.empty( len(inFile), dtype='object')
		errors_pim = np.empty( len(inFile), dtype='object')
			
		hName = key
		k = 0	
		for i in comp_files:
			file = inFile[i]	
			
			hist_pip[i] = file[key]
			hist_pim[i] = file[hName.replace("pip", "pim")]
			xEdges = hist_pip[i].axis(0).edges()
			yEdges = hist_pip[i].axis(1).edges()
			values_pip[i] = hist_pip[i].values()
			values_pip[i][values_pip[i] == 0] = np.nan
			values_pip[i] = np.reshape( values_pip[i], (len(xEdges) - 1, len(yEdges) - 1) )

			values_pim[i] = hist_pim[i].values()
			values_pim[i][values_pim[i] == 0] = np.nan
			values_pim[i] = np.reshape( values_pim[i], (len(xEdges) - 1, len(yEdges) - 1) )

			for j in range(2):
				axs[k,j].xaxis.set_minor_locator(ticker.AutoMinorLocator())
				axs[k,j].yaxis.set_minor_locator(ticker.AutoMinorLocator())
				axs[k,j].tick_params(which='both', width=2)
				axs[k,j].tick_params(which='major', length=7)
				axs[k,j].tick_params(which='minor', length=4)
				setTitle2D( axs[k,j], key)

			axs[k, 0].pcolormesh(xEdges, yEdges, values_pip[i].T, shading='auto')
			axs[k, 1].pcolormesh(xEdges, yEdges, values_pim[i].T,  shading='auto')
			
			for ax in axs.flat:
				ax.label_outer()

			k = k+1
		#fig.legend(fontsize=16)

		fig.savefig(f'comp_plots/{hName}.pdf')
			
