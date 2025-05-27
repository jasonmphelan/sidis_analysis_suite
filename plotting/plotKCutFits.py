import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys

fileName = sys.argv[1] #"../histograms/analysis_note/k2pi_fits.root"
outBase = sys.argv[2]
inFile = uproot.open(fileName)
inFile_R = ROOT.TFile.Open(fileName, "READ")
key = inFile.keys()

def eval_func(params, xVal):
	return params[0]*math.exp(-.5*((xVal - params[1])/params[2])*((xVal-params[1])/params[2]))

chargeSt = ['p','m']
titSt = [r"($(e, e' h^+)$",r"$(e, e' h^-)$"]
for charge in chargeSt:
	for p in range(4):
		if p < 3:
			continue
		for x in range( 14 ):
			
			for q in range( 12 ):
				chIdx = ""
				if charge == 'pip':
					chIdx = 0
				else:
					chIdx = 1
				fig, axs = plt.subplots(7, 4, figsize=(28,32), layout='constrained')
				fig.suptitle(titSt[chIdx] + r"${0:.1f} < Q^2 < {1:.1f}, {2:.2f} < x_B < {3:.2f}$".format(2 + .5*q, 2 + .5*(q+1), .1 + .04*x, .1+0.04*(x+1)), fontsize=14)
				makePlot = False
				for z in range( 28 ):
					reg = z%7
					ch = math.floor(z/7)
					axs[reg, ch].set_title(r"${0:.3f} < z < {1:.3f}$".format( .3 + .025*z, .3 + .025*(z+1)), fontsize=14 )
					axs[reg, ch].set_ylabel(r"Counts [a.u.]", fontsize=14)
					axs[reg, ch].set_xlabel(r"$M^{2}_{RICH}$ [GeV$^2$]", fontsize=14)
					if f"hBeta_rich_pi{charge}_Q2_{q+1}_xB_{x+1}_Z_{z+1}_p_{p+1};1" not in key:
						continue
					if f"f_pi{charge}_{q+1}_{x+1}_{z+1}_{p+1};1" not in key:
						continue
					makePlot = True
					hist = inFile[f"hBeta_rich_pi{charge}_Q2_{q+1}_xB_{x+1}_Z_{z+1}_p_{p+1}"]
					fFit = inFile_R.Get(f"f_pi{charge}_{q+1}_{x+1}_{z+1}_{p+1}")
					
					fit = [fFit.GetParameter(i) for i in range(3)]
					values = hist.values()
					values[values == 0] = np.nan
					errors = hist.errors()
					errors[errors==0] = np.nan

					xEdges = hist.axis().edges()
					binCenters = (xEdges[:-1]+xEdges[1:])/2

					xVals = np.linspace(0, .4, 1000)
					f_vals = [eval_func(fit, xVal) for xVal in xVals]

					axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
					axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

					axs[reg, ch].tick_params(which='both', width=2)
					axs[reg, ch].tick_params(which='major', length=7)
					axs[reg, ch].tick_params(which='minor', length=4)
					
					#not showing errorbars for visual clarity
					#axs[reg, ch].errorbar(binCenters, values, errors, capsize=2, fmt='b', ecolor='b')
					axs[reg, ch].plot(binCenters, values, c='b' ,drawstyle='steps-mid')

					fitProb = fFit.GetProb()
					if fit[1] < -0.05 or fit[1] > 0.075 or fitProb >= 1 or fitProb <= 0 or np.nansum(values) <= 20: 
						axs[reg, ch].axvline( x = 0.1, ymin = 0, ymax = 1, c = 'r', linestyle='-', linewidth=3)
					else:
						axs[reg, ch].plot( xVals, f_vals, c='r', linestyle='-', linewidth=3)

					#if ch != 0:
					#	axs[reg, ch].set_ylabel("")
					#	axs[reg, ch].set_yticks(color='w')

				#plt.show()
				if makePlot:
					#fig.savefig(f"kaon_fits/fits_k2pi_{charge}_p_{p}_xB_{x}_q2_{q}.pdf")
					fig.savefig(outBase+f"/fits_k2pi_{charge}_p_{p}_xB_{x}_q2_{q}.pdf")
				plt.close(fig)
