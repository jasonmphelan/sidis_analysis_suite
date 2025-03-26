import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math

fileName_0 = ""
fileName = "test.root"
inFile = uproot.open(fileName)
inFile_R = ROOT.TFile.Open(fileName, "READ")
key = inFile.keys()

def eval_func(params, xVal):
	return params[0] + params[1]*xVal + params[2]*xVal*xVal

chargeSt = ['p','m']
titSt = [r"($(e, e' h^+)$",r"$(e, e' h^-)$"]
for charge in chargeSt:
	for x in range( 14 ):
		
		chIdx = ""
		if charge == 'pip':
			chIdx = 0
		else:
			chIdx = 1
		fig, axs = plt.subplots(3, 4, figsize=(28,32), layout='constrained')
		#fig.suptitle(titSt[chIdx]+r"${2:.2f} < x_B < {3:.2f}$".format( .1 + .04*x, .1+0.04*(x+1)), fontsize=14)
		makePlot = False
		for q in range( 12 ):
			reg = q%3
			ch = math.floor(q/3)
			axs[reg, ch].set_title(r"${0:.1f} < Q^2 < {1:.1f}$".format( .3 + .025*q, .3 + .025*(q+1)), fontsize=14 )
			axs[reg, ch].set_ylabel(r"Correction", fontsize=14)
			axs[reg, ch].set_xlabel(r"$x$", fontsize=14)
			if f"hMcCorrection_{charge}_1d_{q}_{x};1" not in key:
				print(f"kaonCorr_{charge}_1d_{q+1}_{x+1};1")
				continue
			if f"fitPi{charge}_{q}_{x};1" not in key:
				print(f"fitPi{charge}_{q+1}_{x+1};1")
				continue
			makePlot = True
			hist = inFile[f"hMcCorrection_{charge}_1d_{q}_{x}_smooth"]
			hist_0 = inFile[f"hMcCorrection_{charge}_1d_{q}_{x}"]
		
			fFit =inFile_R.Get(f"fitPi{charge}_{q}_{x}" )
			
			fit = [fFit.GetParameter(i) for i in range(3)]
			values = hist.values()
			values[values == 0] = np.nan
			errors = hist.errors()
			errors[errors==0] = np.nan
			
			values_0 = hist_0.values()
			values_0[values_0 == 0] = np.nan
			errors_0 = hist_0.errors()
			errors_0[errors_0==0] = np.nan

			xEdges = hist.axis().edges()
			binCenters = (xEdges[:-1]+xEdges[1:])/2

			xVals = np.linspace(0.3, .8, 1000)
			f_vals = [eval_func(fit, xVal) for xVal in xVals]

			axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

			axs[reg, ch].tick_params(which='both', width=2)
			axs[reg, ch].tick_params(which='major', length=7)
			axs[reg, ch].tick_params(which='minor', length=4)
			
			#not showing errorbars for visual clarity
			#axs[reg, ch].errorbar(binCenters, values, errors, capsize=2, fmt='b', ecolor='b')
			axs[reg, ch].plot(binCenters, values, c='b' ,drawstyle='steps-mid', label = 'smoothed')
			axs[reg, ch].plot(binCenters, values_0, c='r' ,drawstyle='steps-mid', label = 'Raw Corrections')

			axs[reg, ch].plot( xVals, f_vals, c='k', linestyle='-', linewidth=3)
			axs[reg, ch].set_ylim( 0, 8 )
			axs[reg,ch].legend()

		#plt.show()
		if makePlot:
			fig.savefig(f"mc_fits/smooth_{charge}_xB_{x}.pdf")
		plt.close(fig)
