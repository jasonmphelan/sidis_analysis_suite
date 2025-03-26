import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math

fileName = "../data/correctionFiles/corrections_pi2k_fit.root.root"
inFile = uproot.open(fileName)
inFile_R = ROOT.TFile.Open(fileName, "READ")
key = inFile.keys()

def eval_func(params, xVal):
	#return params[0]*math.exp(-.5*((xVal - params[1])/params[2])*((xVal-params[1])/params[2]))
	#return params[0] + params[1]*xVal + params[2]*xVal**2 + params[3]*xVal**3 
	return params[0]/(1+math.exp(-params[2]*(x-params[3])))+params[1]

chargeSt = ['p','m']
titSt = [r"($(e, e' h^+)$",r"$(e, e' h^-)$"]
for charge in chargeSt:
	for p in range(4):
		if p < 3:
			continue
		for x in range( 14 ):
			chIdx = ""
			if charge == 'pip':
				chIdx = 0
			else:
				chIdx = 1
			fig, axs = plt.subplots(3, 4, figsize=(28,32), layout='constrained')
			fig.suptitle(titSt[chIdx] + r", {0} < x_B < {1}$".format( .1 + .04*x, .1+0.04*(x+1)), fontsize=14)
			makePlot = False
			
			for q in range( 12 ):
				reg = q%3
				ch = math.floor(q/3)
				axs[reg, ch].set_title(r"${0} < Q^2 < {1}$".format( 2 + .5*q, 2 + .5*(q+1)) )
				axs[reg, ch].set_ylabel(r"w", fontsize=14)
				axs[reg, ch].set_xlabel(r"$z$", fontsize=14)
				if f"kaonCorr_{charge}_1d_{p}_{q}_{x};1" not in key:
					continue
				if f"fitPi{charge}_{p}_{q}_{x};1" not in key:
					continue
				makePlot = True
				hist = inFile[f"kaonCorr_{charge}_1d_{p}_{q}_{x}"]
				fFit = inFile_R.Get(f"fitPi{charge}_{p}_{q}_{x}")
				

				fit = [fFit.GetParameter(i) for i in range(4)]
				values = hist.values()
				values[values == 0] = np.nan
				errors = hist.errors()
				errors[errors==0] = np.nan

				xEdges = hist.axis().edges()
				binCenters = (xEdges[:-1]+xEdges[1:])/2

				#xVals = np.linspace(0, .4, 1000)
				xVals = np.linspace(0, 1, 1000)
				f_vals = [eval_func(fit, xVal) for xVal in xVals]

				axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
				axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

				axs[reg, ch].tick_params(which='both', width=2)
				axs[reg, ch].tick_params(which='major', length=7)
				axs[reg, ch].tick_params(which='minor', length=4)
				
				#not showing errorbars for visual clarity
				axs[reg, ch].errorbar(binCenters, values, errors, capsize=2, fmt='b', ecolor='b', linestyle='')
				#axs[reg, ch].plot(binCenters, values, c='b' ,drawstyle='steps-mid')

				axs[reg, ch].plot( xVals, f_vals, c='r', linestyle='-', linewidth=3)
				axs[reg, ch].set_ylim( [ 0, 1.2 ] )				
				axs[reg, ch].set_xlim( [ 0.3, 1 ] )				

				#if ch != 0:
				#	axs[reg, ch].set_ylabel("")
				#	axs[reg, ch].set_yticks(color='w')

			#plt.show()
			if makePlot:
				fig.savefig(f"kaon_corrections_fits/fits_pi2k_{charge}_p_{p}_xB_{x}.pdf")
			plt.close(fig)
