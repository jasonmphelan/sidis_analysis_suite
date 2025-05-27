import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys

in_name = sys.argv[1]
out_dir = sys.argv[2]
inFile = uproot.open(in_name)
 
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
		plt.xlabel(r"$\beta$", fontsize=18)
	if "Vz_e" in key:
		plt.xlabel(r"$v_z^e$ [cm]", fontsize=18)
	if "Vz_pi" in key:
		plt.xlabel(r"$v_z^{\pi} - v_z^e$ [cm]", fontsize=18)
	if "Z_pi" in key:	
		plt.xlabel(r"$z$", fontsize=18)
	if "Y_pi" in key:
		plt.xlabel(r"$y$", fontsize=18)
	if "W_pi" in key:
		plt.xlabel(r"$W$ [GeV]", fontsize=18)
	if "Theta_e" in key:
		plt.xlabel(r"$\theta_e$ [deg.]", fontsize=18)
	if "Theta_pi" in key:
		plt.xlabel(r"$\theta_{\pi} [deg.]$", fontsize=18)
	if "Q2" in key:
		plt.xlabel(r"$Q^2$ [GeV$^2$]", fontsize=18)
	if "P_pi" in key:
		plt.xlabel(r"$p_{\pi}$ [GeV]", fontsize=18)
	if "P_e" in key:
		plt.xlabel(r"$p_{e}$ [GeV]", fontsize=18)
	if "Mx" in key:
		plt.xlabel(r"$M_X$ [GeV]", fontsize=18)
	if "Omega" in key:
		plt.xlabel(r"$\omega$ [GeV]", fontsize=18)


#inFile = uproot.open("../histograms/analysis_note/kinematics_data_10.4.root")
for key in inFile.keys():
	if "pim" in key:
		continue
	hist = inFile[key]

	if isinstance( hist, uproot.models.TH.Model_TH1F_v3): 
		fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		pip = inFile[key]
		hName = key
		pim = inFile[hName.replace("pip", "pim")]

		values_pip = pip.values()		
		errors_pip = pip.errors()
		values_pim = pim.values()
		errors_pim = pim.errors()

		values_pip[values_pip == 0] = np.nan
		values_pim[values_pim == 0] = np.nan

		ratio = values_pip/values_pim
		ratio_err = ratio*np.sqrt( (errors_pip/values_pip)**2 + (errors_pim/values_pim)**2)


		values_pip = values_pip/np.nansum(pip.values())
		errors_pip = errors_pip/np.nansum(pip.values())
		values_pim = values_pim/np.nansum(pim.values())
		errors_pim = errors_pim/np.nansum(pim.values())

		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2

		ax[0].errorbar(binCenters, values_pip, errors_pip, capsize=2, fmt='r', ecolor='r',drawstyle='steps-mid', label=r"$(e, e'\pi^+)$")
		#axs[0].set_title(rf"$d(e, e'\pi^+)$")

		#drawCut1D( axs[0], key )		
		setTitle1D( fig, key )

		ax[0].errorbar(binCenters, values_pim, errors_pim, capsize=2, fmt='b', ecolor='b',drawstyle='steps-mid', label=r"$(e, e'\pi^-)$")
		
		ax[1].errorbar(binCenters, ratio, ratio_err, capsize=2, fmt='k o', ecolor='k')
		
		ax[1].axhline(y=1, linestyle='--', color='grey')	
		
		plt.margins(y=0)
		plt.margins(x=0)
		
		for i in range(2):
			ax[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			ax[i].tick_params(which='both', width=2, labelsize=12)
			ax[i].tick_params(which='major', length=7, labelsize=12)
			ax[i].tick_params(which='minor', length=4, labelsize=12)
		
		ax[0].set_ylabel("Counts [a.u.]", fontsize=18)
		ax[1].set_ylabel(r"$Y(e, e'\pi^+)/Y(e, e'\pi^-)$", fontsize=12)
		#setTitle1D( ax[1], key )

		pdfName = hName.replace("_pip", "")
		fig.legend(fontsize=16)
		fig.savefig(f"{out_dir}/{pdfName}.pdf")
		plt.close(fig)
