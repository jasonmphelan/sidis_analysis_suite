import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
plt.rcParams["font.family"] = "sans-serif"

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

in_names = int(sys.argv[1])
out_dir = int(sys.argv[2])
inFile = uproot.open(in_names[0])#"../histograms/analysis_note/kinematics_data_10.4.root")
inFile_MC = uproot.open(in_names[1])#"../histograms/analysis_note/kinematics_gemc_10.4.root")
for key in inFile.keys():
	if "pim" in key:
		continue
	
	pip = inFile[key]
	if isinstance( pip, uproot.models.TH.Model_TH1F_v3): 
		fig, ax = plt.subplots(2, 2, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
	        
		pip_MC = inFile_MC[key]
		hName = key
		pim = inFile[hName.replace("pip", "pim")]
		pim_MC = inFile_MC[hName.replace("pip", "pim")]

		values_pip = pip.values()		
		errors_pip = pip.errors()
		values_pip_MC = pip_MC.values()		
		errors_pip_MC = pip_MC.errors()
		
		values_pim = pim.values()
		errors_pim = pim.errors()
		values_pim_MC = pim_MC.values()		
		errors_pim_MC = pim_MC.errors()

		values_pip[values_pip == 0] = np.nan
		values_pim[values_pim == 0] = np.nan
		values_pip_MC[values_pip_MC == 0] = np.nan
		values_pim_MC[values_pim_MC == 0] = np.nan
		
		values_pip = values_pip/np.nansum(pip.values())
		errors_pip = errors_pip/np.nansum(pip.values())
		values_pim = values_pim/np.nansum(pim.values())
		errors_pim = errors_pim/np.nansum(pim.values())
		
		values_pip_MC = values_pip_MC/np.nansum(pip_MC.values())
		errors_pip_MC = errors_pip_MC/np.nansum(pip_MC.values())
		values_pim_MC = values_pim_MC/np.nansum(pim_MC.values())
		errors_pim_MC = errors_pim_MC/np.nansum(pim_MC.values())

		ratio_pip = values_pip/values_pip_MC
		ratio_pim = values_pip_MC/values_pim_MC
		
		ratio_err_pip = ratio_pip*np.sqrt( (errors_pip_MC/values_pip_MC)**2 + (errors_pip/values_pip)**2)
		ratio_err_pim = ratio_pim*np.sqrt( (errors_pim_MC/values_pim_MC)**2 + (errors_pim/values_pim)**2)



		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2

		ax[0, 0].errorbar(binCenters, values_pip, errors_pip, capsize=2, fmt='r', ecolor='r',drawstyle='steps-mid', label="Data")
		ax[0, 0].errorbar(binCenters, values_pip_MC, errors_pip, capsize=2, fmt='b o', ecolor='b', label="GEMC")
		ax[0,0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
		ax[0,0].set_ylabel('Counts [a.u.]', fontsize=16)
		ax[1, 0].errorbar(binCenters, ratio_pip, ratio_err_pip, capsize=2, fmt='k o', ecolor='k')
		setTitle1D( ax[1,0], key )	
		ax[1, 0].set_ylabel(r'$Y_{Data}/Y_{GEMC}$', fontsize=16)		

		ax[0, 1].errorbar(binCenters, values_pim, errors_pim, capsize=2, fmt='r', ecolor='r',drawstyle='steps-mid')
		ax[0, 1].errorbar(binCenters, values_pim_MC, errors_pim, capsize=2, fmt='b o', ecolor='b')
		ax[0,1].set_title(rf"$d(e, e'\pi^-)$", fontsize=16)
		#setTitle1D( fig, key )
		ax[1, 1].errorbar(binCenters, ratio_pim, ratio_err_pim, capsize=2, fmt='k o', ecolor='k')
		setTitle1D( ax[1,1], key )
		
		

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
		fig.savefig(f"{out_dir}/{pdfName}.pdf")
		plt.close(fig)
