import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys

fileName = sys.argv[1] #"../histograms/analysis_note/k2pi_fits.root"
fileName_2 = sys.argv[2]
outBase = sys.argv[3]
inFile = uproot.open(fileName)
inFile_R = ROOT.TFile.Open(fileName_2, "READ")
inFile_F = uproot.open(fileName_2)
key = inFile.keys()
key_f = inFile_F.keys()

def eval_func(params, xVal):
	return params[0]*math.exp(-.5*((xVal - params[1])/params[2])*((xVal-params[1])/params[2]))

chargeSt = ['p','m']
titSt = [r"($(e, e' h^+)$",r"$(e, e' h^-)$"]
p_bins = [1.25, 2, 3, 4, 5]
for charge in chargeSt:
	for p in range(4):
		if p < 2:
			continue
		for x in range( 7 ):
			
			for q in range( 12 ):
				chIdx = ""
				if charge == 'pip':
					chIdx = 0
				else:
					chIdx = 1
				fig, axs = plt.subplots(3, 4, figsize=(28,32), layout='constrained')
				fig.suptitle(titSt[chIdx] + r"${0:.1f} < Q^2 < {1:.1f}, {2:.2f} < x_B < {3:.2f}, {4} < p < {5}$".format(2 + .5*q, 2 + .5*(q+1), .1 + .08*x, .1+0.08*(x+1), p_bins[p], p_bins[p+1]), fontsize=22)
				makePlot = False
				for z in range( 12 ):
					ch = z%4
					reg = math.floor(z/4)
					axs[reg, ch].set_title(r"${0:.2f} < z < {1:.2f}$".format( .3 + .05*z, .3 + .05*(z+1)), fontsize=22 )
					axs[reg, ch].set_ylabel(r"Counts [a.u.]", fontsize=22)
					axs[reg, ch].set_xlabel(r"$M^{2}_{RICH}$ [GeV$^2$]", fontsize=22)
					axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
					axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())
					axs[reg, ch].tick_params(which='both', width=2, labelsize=20)
					axs[reg, ch].tick_params(which='major', length=7)
					axs[reg, ch].tick_params(which='minor', length=4)
				
					if f"hBeta_rich_pi{charge}_Q2_{q+1}_xB_{x+1}_Z_{z+1}_p_{p+1};1" not in key:
						axs[reg, ch].text( 0.45, .5, "Empty", fontsize=30)
						continue

					makePlot = True
					hist = inFile[f"hBeta_rich_pi{charge}_Q2_{q+1}_xB_{x+1}_Z_{z+1}_p_{p+1}"]
					xVals = []
					f_vals = []
					
					values = hist.values()
					
					values[values == 0] = np.nan
					errors = hist.errors()
					errors[errors==0] = np.nan

					xEdges = hist.axis().edges()
					binCenters = (xEdges[:-1]+xEdges[1:])/2

					#xVals = np.linspace(-.1, .4, 1000)
					#f_vals = [eval_func(fit, xVal) for xVal in xVals]

					axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
					axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

					axs[reg, ch].tick_params(which='both', width=2)
					axs[reg, ch].tick_params(which='major', length=7)
					axs[reg, ch].tick_params(which='minor', length=4)
					
					#not showing errorbars for visual clarity
					#axs[reg, ch].errorbar(binCenters, values, errors, capsize=2, fmt='b', ecolor='b')
					axs[reg, ch].plot(binCenters, values, c='b' ,drawstyle='steps-mid')

					#axs[reg, ch].text( 0.1, 0.75*np.max(values), "Counts = {0}".format(np.sum(values)))

					if f"f_pi{charge}_{q+1}_{x+1}_{z+1}_{p+1};1" in key_f:
						#continue
						fFit = inFile_R.Get(f"f_pi{charge}_{q+1}_{x+1}_{z+1}_{p+1}")
					
						fit = [fFit.GetParameter(i) for i in range(3)]
						xVals = np.linspace(-.1, .4, 1000)
						f_vals = [eval_func(fit, xVal) for xVal in xVals]

						fitProb = fFit.GetProb()
						#axs[reg, ch].text( 0.2, 0.75*np.max(values), "P-Value = {0:.2f}".format(fFit.GetProb()), fontsize=24)
						if fit[1] < -0.05 or fit[1] > 0.075 or fitProb >= 1 or fitProb <= 0 or np.nansum(values) <= 150 or fit[1] + 2.5*fit[2]>0.2: 
							axs[reg, ch].axvline( x = 0.1, ymin = 0, ymax = 1, c = 'r', linestyle='-', linewidth=3)
							axs[reg, ch].text( 0.11, 0.85*np.nanmax(values), "Counts = {0}".format(np.nansum(values)), fontsize=24)
						else:
							#fitProb = fFit.GetProb()
							
							axs[reg, ch].text( 0.2, 0.8*np.nanmax(values), "P-Value = {0:.2f}".format(fFit.GetProb()), fontsize=24)
							axs[reg, ch].plot( xVals, f_vals, c='r', linestyle='-', linewidth=3)
							axs[reg, ch].axvline( x = fit[1] + 2.5*fit[2], ymin = 0, ymax = 1, c = 'k', linestyle='--', linewidth=3)
							axs[reg, ch].text( 0.2, 0.85*np.nanmax(values), "Counts = {0}".format(np.nansum(values)), fontsize=24)
					else:
						axs[reg, ch].axvline( x = 0.1, ymin = 0, ymax = 1, c = 'r', linestyle='-', linewidth=3)
						axs[reg, ch].text( 0.11, 0.85*np.nanmax(values), "Counts = {0}".format(np.nansum(values)), fontsize=24)

					#if ch != 0:
					#	axs[reg, ch].set_ylabel("")
					#	axs[reg, ch].set_yticks(color='w')

				#plt.show()
				if makePlot:
					#fig.savefig(f"kaon_fits/fits_k2pi_{charge}_p_{p}_xB_{x}_q2_{q}.pdf")
					if 'pi2k' in fileName:
						fig.savefig(outBase+f"/fits_pi2k_{charge}_p_{p}_xB_{x}_q2_{q}.pdf")
					else:
						fig.savefig(outBase+f"/fits_k2pi_{charge}_p_{p}_xB_{x}_q2_{q}.pdf")
				plt.close(fig)
