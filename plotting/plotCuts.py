import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import sys 

in_name = (sys.argv[1])
out_dir = (sys.argv[2])

def setTitle2D( ax, key):
	if "Beta" in key or "Min" in key or "Pos" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=18)
		ax.set_ylabel(r"$\beta$", fontsize=18)
	else:
		if "WV" in key:
			ax.set_xlabel(r"V$_{PCAL}$ [cm]", fontsize=18)
			ax.set_ylabel(r"W$_{PCAL}$ [cm]", fontsize=18)
		if "Edep" in key:
			ax.set_xlabel(r"$E_{PCAL}$ [GeV]", fontsize=18)
			ax.set_ylabel(r"$E_{CIN} + E_{COUT}$ [GeV]", fontsize=18)
		if "SF_corr" in key:
			ax.set_xlabel(r"$\Delta E_{dep}(PCAL)/p_{e}$", fontsize=18)
			ax.set_ylabel(r"$(E_{PCAL} + E_{CIN} + E_{COUT})/p_e$", fontsize=18)
		if "SF_sec" in key:
			ax.set_xlabel(r"$p_{e}$ [GeV]", fontsize=18)
			ax.set_ylabel(r"$\Delta E_{CIN}/p_e$", fontsize=18)
		if "Chi2" in key:
			ax.set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=18)
			ax.set_ylabel(r"$\chi^2_{PID}$", fontsize=18)
		if "Q2_omega" in key:
			ax.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=18)
			ax.set_ylabel(r"$\omega$ [GeV]", fontsize=18)
		if "Q2_W" in key:
			ax.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=18)
			ax.set_ylabel(r"$W$ [GeV]", fontsize=18)
		if "Q2_Z" in key:
			ax.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=18)
			ax.set_ylabel(r"$z$ [GeV]", fontsize=18)
def setTitle1D( ax, key):
	ax.set_ylabel("Counts [a.u.]", fontsize=18)
	if "Beta" in key or "Min" in key or "Pos" in key:
		ax.set_xlabel(r"$\beta$", fontsize=18)
		ax.set_yscale("log")
	else:
		if "Vz_e" in key:
			ax.set_xlabel(r"$v_z^e$ [cm]", fontsize=18)
		if "Vz_pi" in key:
			ax.set_xlabel(r"$v_z^{\pi} - v_z^e$ [cm]", fontsize=18)
		if "Z_pi" in key:	
			ax.set_xlabel(r"$z$", fontsize=18)
		if "Y_pi" in key:
			ax.set_xlabel(r"$y$", fontsize=18)
		if "W_pi" in key:
			ax.set_xlabel(r"$W$ [GeV]", fontsize=18)
		if "Theta_e" in key:
			ax.set_xlabel(r"$\theta_e$ [deg.]", fontsize=18)
		if "Theta_pi" in key:
			ax.set_xlabel(r"$\theta_{\pi} [deg.]$", fontsize=18)
		if "Q2" in key:
			ax.set_xlabel(r"$Q^2$ [GeV$^2$]", fontsize=18)
		if "P_pi" in key:
			ax.set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=18)
		if "P_e" in key:
			ax.set_xlabel(r"$p_{e}$ [GeV]", fontsize=18)
		if "Mx" in key:
			ax.set_xlabel(r"$M_X$ [GeV]", fontsize=18)


def drawCut2D( ax, key ):
	
	if "Beta" in key or "Min" in key or "Pos" in key:
		pass
	else:
		if "WV" in key:
			ax.axvline(x=19, ymin=19/200, ymax=1, c="red", linewidth=2)
			ax.axhline(y=19, xmin=19/200, xmax=1, c="red", linewidth=2)
		if "Edep" in key:
			ax.axvline(x=0.07, c="red", linewidth=2)
		if "SF_corr" in key:
			ax.plot([0,.2],[.2,0], 'r-', linewidth=2)
		if "Chi2" in key:
			C = .93
			if "pip" in key:
				C = .88
			p1=np.arange(2.44, 4.6, .01)
			p2=np.arange(4.6, 9, .01)
			ax.axhline(y = -3*C, xmin=0, xmax=1, c='red', linewidth=2)
			ax.axhline(y = 3*C, xmin=0, xmax=2.44/9, c='red', linewidth=2)
			ax.plot( p1, C*( 0.00869 + 14.98587*np.exp(-p1/1.18236) + 1.81751*np.exp(-p1/4.86394) ), 'r-', linewidth=2)
			ax.plot( p2, C*(-1.14099 + 24.14992*np.exp(-p2/1.3655)+2.66876*np.exp(-p2/6.80552) ), 'r-', linewidth=2)
		if "Q2_W" in key:
			ax.axvline(x=2, ymin = 1/3, ymax = 1, c='r', linewidth=2)
			ax.axhline(y=2.5, xmin = 2/10, xmax=1, c='r', linewidth=2)

def drawCut1D( ax, key ):
	if "Beta" in key or "Min" in key or "Pos" in key:
		pass
	else:
		if "Epcal" in key or "EPcal" in key:
			ax.axvline(x=.07, c="red", linewidth=2) 
		if "Vz_e" in key:
			ax.axvline(x=-7., c='r', linewidth=2)
			ax.axvline(x=2, c='r', linewidth=2)
		if "Vz_pip" in key:
			ax.axvline(x=-7, c='r', linewidth=2)
			ax.axvline(x= 4, c='r', linewidth=2)
		if "Vz_pim" in key:
			ax.axvline(x=-5, c='r', linewidth=2)
			ax.axvline(x= 4, c='r', linewidth=2)
		if "Z_pi" in key:	
			ax.axvline(x=.3, c='r', linewidth=2)
		if "Y_pi" in key:
			ax.axvline(x=.75, c='r', linewidth=2)
		if "W_pi" in key:
			ax.axvline(x=2.5, c='r', linewidth=2)
		if "Theta" in key:
			ax.axvline(x=5, c='r', linewidth=2)
			ax.axvline(x=35, c='r', linewidth=2)
		if "Q2" in key:
			ax.axvline(x=2, c='r', linewidth=2)
		if "P_pi" in key:
			ax.axvline(x=1.25, c='r', linewidth=2)
			ax.axvline(x=5, c='r', linewidth=2)
		if "P_e" in key:
			ax.axvline(x=3, c='r', linewidth=2)
		if "Mx" in key:
			ax.axvline(x=1.7, c='r', linewidth=2)


#inFile = uproot.open("../histograms/analysis_note/kinematic_cuts_10.2.root")
#inFile = uproot.open("../histograms/analysis_note/detector_plots_10.4.root")
inFile = uproot.open(in_name)
for key in inFile.keys():
	if "pim" in key:
		continue
	hist = inFile[key]
	if isinstance( hist, uproot.models.TH.Model_TH2F_v4): 
		if "reg" in key:
			continue
		fig, axs = plt.subplots(1, 2, figsize=(18,6), layout='constrained')
		pip = inFile[key]
		hName = key
		pim = inFile[hName.replace("pip", "pim")]

		values_pip = pip.values()
		values_pip[values_pip == 0] = np.nan
		values_pim = pim.values()
		values_pim[values_pim == 0] = np.nan

		xEdges = pip.axis(0).edges()
		yEdges = pim.axis(1).edges()
		
		axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
	
		axs[0].tick_params(which='both', width=2)
		axs[0].tick_params(which='major', length=7)
		axs[0].tick_params(which='minor', length=4)
		axs[1].tick_params(which='both', width=2)
		axs[1].tick_params(which='major', length=7)
		axs[1].tick_params(which='minor', length=4)
		
		axs[0].pcolormesh(xEdges, yEdges, values_pip.T, cmap='viridis', shading='auto', norm="log")
		axs[0].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)

		drawCut2D( axs[0], key )		
		setTitle2D( axs[0], key )

		axs[1].pcolormesh(xEdges, yEdges, values_pim.T, cmap='viridis', shading='auto', norm="log")
		axs[1].set_title(rf"$(e, e'\pi^-)$", fontsize=16)
		
		drawCut2D( axs[1], hName.replace("pip","pim") )		
		setTitle2D( axs[1], key )

		pdfName = hName.replace(";1", "")
		pdfName = pdfName.replace("_pip", "")
		fig.savefig(f"{out_dir}/{pdfName}.pdf")
		plt.close(fig)

	if isinstance( hist, uproot.models.TH.Model_TH1F_v3): 
		fig, axs = plt.subplots(1, 2, figsize=(18,6), layout='constrained')
		pip = inFile[key]
		hName = key
		pim = inFile[hName.replace("pip", "pim")]

		values_pip = pip.values()
		errors_pip = pip.errors()
		values_pim = pim.values()
		errors_pim = pim.errors()

		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2
		
		axs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
	
		axs[0].tick_params(which='both', width=2, labelsize=12)
		axs[0].tick_params(which='major', length=7, labelsize=12)
		axs[0].tick_params(which='minor', length=4, labelsize=12)
		axs[1].tick_params(which='both', width=2, labelsize=12)
		axs[1].tick_params(which='major', length=7, labelsize=12)
		axs[1].tick_params(which='minor', length=4, labelsize=12)

		#axs[0].plot(binCenters, values_pip, 'k', drawstyle='steps-mid')
		axs[0].errorbar(binCenters, values_pip, errors_pip, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid')
		axs[0].set_title(rf"$d(e, e'\pi^+)$",fontsize=16)

		drawCut1D( axs[0], key )		
		setTitle1D( axs[0], key )

		#axs[1].plot(binCenters, values_pim, 'k', drawstyle='steps-mid')
		axs[1].errorbar(binCenters, values_pim, errors_pim, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid')
		axs[1].set_title(rf"$d(e, e'\pi^-)$",fontsize=16)

		pimName = hName.replace("_pip", "_pim")	
	
		drawCut1D( axs[1], pimName )
		setTitle1D( axs[1], key )

		pdfName = hName.replace("_pip", "")
		pdfName = pdfName.replace(";1", "")
		fig.savefig(f"{out_dir}/{pdfName}.pdf")
		plt.close(fig)
