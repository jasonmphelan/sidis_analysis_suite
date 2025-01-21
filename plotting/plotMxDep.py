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
def setTitle1D( ax, key):
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
		ax.set_xlabel(r"$M_X$ [GeV]", fontsize=16)
	if "Omega" in key:
		ax.set_xlabel(r"$\omega$ [GeV]", fontsize=16)


inFile = uproot.open("../histograms/analysis_note/Mx_dep.root")

for key in inFile.keys():
	if "pim" in key:
		continue
	if "Mx" in key:
		continue

	pip = inFile[key]
	fig, ax = plt.subplots(2, 2, figsize=(12,6), height_ratios=[2,4], layout='constrained')
	
	hName = key
	pim = inFile[hName.replace("pip", "pim")]

	values_pip = pip.values()		
	errors_pip = pip.errors()
	
	values_pim = pim.values()
	errors_pim = pim.errors()

	values_pip[values_pip == 0] = np.nan
	values_pim[values_pim == 0] = np.nan
	
	xEdges = pip.axis().edges()
	binCenters = (xEdges[:-1]+xEdges[1:])/2

	ax[0, 0].errorbar(binCenters, values_pip, errors_pip, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid', label="Data")
	ax[0,0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
	setTitle1D( ax[0,0], key )	
	ax[1,1].set_ylabel('Counts [a.u.]', fontsize=16)
	ax[0,0].set_ylabel('Counts [a.u.]', fontsize=16)

	ax[0, 1].errorbar(binCenters, values_pim, errors_pim, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid')
	ax[0,1].set_title(rf"$d(e, e'\pi^-)$", fontsize=16)
	ax[0,1].set_ylabel('Counts [a.u.]', fontsize=16)
	setTitle1D( ax[0,1], key )
	ax[0,1].set_yscale("log")
	ax[0,0].set_yscale("log")
	
	keyBase = ''
	if 'Z' in hName:
		keyBase = 'hMx_Zpip_'
	if 'Q2' in hName:
		keyBase = 'hMx_Q2pip_'
	ran = [3,2,1,0]

	for b in ran:
		pip = inFile[keyBase+f'{b};1']
	
		pim = inFile[(keyBase+f'{b};1').replace("pip", "pim")]

		values_pip = pip.values()		
		errors_pip = pip.errors()
	
		values_pim = pim.values()
		errors_pim = pim.errors()
		
		values_pip = values_pip/np.nansum(pip.values())
		errors_pip = errors_pip/np.nansum(pip.values())
		values_pim = values_pim/np.nansum(pim.values())
		errors_pim = errors_pim/np.nansum(pim.values())

		values_pip[values_pip == 0] = np.nan
		values_pim[values_pim == 0] = np.nan
	
		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2
	
		ax[1, 0].errorbar(binCenters, values_pip, errors_pip, capsize=1, fmt=colorList[b], ecolor=colorList[b],drawstyle='steps-mid', label="Data")
		ax[1, 0].set_ylabel('Counts [a.u.]', fontsize=16)
		ax[1, 0].set_xlabel(r'$M_X$ [GeV]', fontsize=16)

		ax[1, 1].errorbar(binCenters, values_pim, errors_pim, capsize=1, fmt=colorList[b], ecolor=colorList[b],drawstyle='steps-mid')
		ax[1, 1].set_ylabel('Counts [a.u.]', fontsize=16)
		ax[1, 1].set_xlabel(r'$M_X$ [GeV]', fontsize=16)
		

	ax[1,0].axvline(1.7, c='k', linestyle='--')
	ax[1,1].axvline(1.7, c='k', linestyle='--')

	plt.margins(y=0)
	plt.margins(x=0)
	
	for j in range(2):
		for i in range(2):
			ax[j,i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			ax[j,i].tick_params(which='both', width=2, labelsize=12)
			ax[j,i].tick_params(which='major', length=7, labelsize=12)
			ax[j,i].tick_params(which='minor', length=4, labelsize=12)
			ax[j,i].margins(y=0)
			ax[j,i].margins(x=0)

	cuts = []
	if "Z" in key:
		cuts = [.3,.4,.5,.6,.8]
	else:
		cuts = [2,3,4,6,8]

	for line in ran:
		ax[0,0].axvspan( cuts[line], cuts[line+1], alpha = .5, color = colorList[line] )
		ax[0,1].axvspan( cuts[line], cuts[line+1], alpha = .5, color = colorList[line] )


	
	#ax[0].set_ylabel("Counts [a.u.]", fontsize=18)
	#ax[1].set_ylabel(r"$Y(e, e'\pi^+)/Y(e, e'\pi^-)$", fontsize=12)
	#setTitle1D( ax[1], key )

	pdfName = hName.replace("_pip", "")
	pdfName = pdfName.replace(";1", "")
	fig.savefig(f"mx_plots/{pdfName}.pdf")
	plt.close(fig)
