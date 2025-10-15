import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys
plt.rcParams["font.family"] = "sans-serif"
colorList = ['red', 'gold', 'blue', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']

accCharge = [ 3.49779e+07 + 2.90652e+07 + 2.60921e+07, 1.16224e+07,1.17782e+07]

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
	if "Eta" in key:
		ax.set_xlabel(r"$\eta$", fontsize=16)
	if "Xb" in key:
		ax.set_xlabel(r"$x_B$", fontsize=16)
	if "Pt_pi" in key:
		ax.set_xlabel(r"$p^{\perp}_{\pi}$ [GeV]", fontsize=16)

inFile_names = sys.argv[1]
out_names = sys.argv[2]

comp_files = {0,1}
labels = [r'Inbending $e^-$', r'Outbending $e^+$ (Hadron Trigger)', r'Outbending $e^+$ (Electron Trigger)']

inFile = [ uproot.open( inFile_names + f'_{i}.root') for i in ['allE', 'out', 'out_etrig'] ]

for key in inFile[0].keys():
	if "pim" in key:
		continue
	
	if isinstance( inFile[0][key], uproot.models.TH.Model_TH1F_v3): 
		fig, ax = plt.subplots(2, 2, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		values_pip_u = np.empty( len(inFile), dtype='object')
		values_pim_u = np.empty( len(inFile), dtype='object')
		errors_pip_u = np.empty( len(inFile), dtype='object')
		errors_pim_u = np.empty( len(inFile), dtype='object')

		values_pip = np.empty( len(inFile), dtype='object')
		values_pim = np.empty( len(inFile), dtype='object')
		errors_pip = np.empty( len(inFile), dtype='object')
		errors_pim = np.empty( len(inFile), dtype='object')
		
		for i in range(len(inFile)):
			file = inFile[i]	
			
			pip = file[key]
			hName = key
			pim = file[hName.replace("pip", "pim")]

			values_pip_u[i] = pip.values()		
			errors_pip_u[i] = pip.errors()
			
			values_pim_u[i] = pim.values()
			errors_pim_u[i] = pim.errors()

			values_pip_u[i][values_pip_u[i] == 0] = np.nan
			values_pim_u[i][values_pim_u[i] == 0] = np.nan
			
			values_pip[i] = values_pip_u[i]/accCharge[i]
			errors_pip[i] = errors_pip_u[i]/accCharge[i]
			values_pim[i] = values_pim_u[i]/accCharge[i]
			errors_pim[i] = errors_pim_u[i]/accCharge[i]

			ratio_list = [1, 2, 2, 2]

			ratio_pip = values_pip[0]/values_pip[i]
			ratio_pim = values_pim[0]/values_pim[i]
			
			
			ratio_err_pip = ratio_pip*np.sqrt( (errors_pip[i]/values_pip[i])**2 + (errors_pip[0]/values_pip[0])**2)
			ratio_err_pim = ratio_pim*np.sqrt( (errors_pim[i]/values_pim[i])**2 + (errors_pim[0]/values_pim[0])**2)

			values_pip[i] =  values_pip_u[i]/np.nansum(values_pip_u[i])
			values_pim[i] =  values_pim_u[i]/np.nansum(values_pim_u[i])
			errors_pip[i] =  errors_pip_u[i]/np.nansum(values_pip_u[i])
			errors_pim[i] =  errors_pim_u[i]/np.nansum(values_pim_u[i])
			
			xEdges = pip.axis().edges()
			binCenters = (xEdges[:-1]+xEdges[1:])/2

			ax[0, 0].errorbar(binCenters, values_pip[i], errors_pip[i], capsize=2, fmt=colorList[i], ecolor=colorList[i],drawstyle='steps-mid', label=labels[i])
			ax[0,0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
			ax[0,0].set_ylabel('Counts [a.u.]', fontsize=16)
			ax[0, 1].errorbar(binCenters, values_pim[i], errors_pim[i], capsize=2, fmt=colorList[i], ecolor=colorList[i],drawstyle='steps-mid')
			ax[0,1].set_title(rf"$d(e, e'\pi^-)$", fontsize=16)
			if i == 1 or i == 2:
				ax[1, 0].errorbar(binCenters, ratio_pip, ratio_err_pip, capsize=2, fmt=colorList[i], ecolor=colorList[i], marker='o', mec='k', linestyle=' ')
				setTitle1D( ax[1,0], key )	
				ax[1, 0].set_ylabel(r'$Y(E)/Y(E = 10.2$ GeV)', fontsize=14)		
				ax[1, 1].errorbar(binCenters, ratio_pim, ratio_err_pim, capsize=2, fmt=colorList[i], ecolor=colorList[i], marker='o', mec='k',linestyle=' ')
				setTitle1D( ax[1,1], key )	
				ax[1, 0].set_ylabel(r'$Y(E)/Y(E = 10.2$ GeV)', fontsize=12)
				ax[1,0].set_ylim(20, 1500)
				ax[1,1].set_ylim(20, 1500)
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
		fig.savefig(f"{out_names}/{pdfName}.pdf")
		plt.close(fig)
			
