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

def eval_func_rho(params, xVal):
	if params[2] == 0:
		return 0
	return params[0]*math.exp(-.5*((xVal - params[1])/params[2])*((xVal-params[1])/params[2]))# + params[3]*math.exp(-.5*((xVal - params[4])/params[5])*((xVal-params[4])/params[5]))
def eval_func_back(params, xVal):
	if params[5] == 0:
		return 0
	return params[3]*math.exp(-.5*((xVal - params[4])/params[5])*((xVal-params[4])/params[5]))
'''
def eval_func_rho(params, xVal):
	#if params[2] == 0:
	#	return 0
	return params[0]/(3.14)*1./((xVal-params[1])**2 + params[0]**2 )# + params[3]*math.exp(-.5*((xVal - params[4])/params[5])*((xVal-params[4])/params[5]))
def eval_func_back(params, xVal):
	#if params[5] == 0:
	#	return 0
	return params[2]/(3.14)*1./((xVal-params[2])**2 + params[2]**2 )
'''

chargeSt = ['p','m']
titSt = [r"($(e, e' h^+)$",r"$(e, e' h^-)$"]
	
for x in range( 14 ):
		
	for q in range( 12 ):
		fig, axs = plt.subplots(3, 3, figsize=(28,32), layout='constrained', sharey='row', sharex='row')
		fig.suptitle( r"${0:.1f} < Q^2 < {1:.1f}, {2:.2f} < x_B < {3:.2f}$".format(2 + .5*q, 2 + .5*(q+1), .1 + .04*x, .1+0.04*(x+1)), fontsize=18)
		makePlot = False
		
		for charge in chargeSt:
			color = 'b'
			pltIdx = 0
			chIdx = ""
			if charge == 'p':
				chIdx = 0
			else:
				chIdx = 1
				pltIdx = 2
				color = 'r'
			
			for z in range( 9 ):
				ch = z%3 
				reg = math.floor(z/3)# + pltIdx
				axs[reg, ch].set_title(r"${0:.3f} < z < {1:.3f}$".format( .3 + .05*z, .3 + .05*(z+1)), fontsize=18 )
				axs[reg, ch].set_ylabel(r"Counts [a.u.]", fontsize=16)
				axs[reg, ch].set_xlabel(r"$M_{X}^{(e,e'\pi^+\pi^-)}$ [GeV]", fontsize=16)
				if f"hMx_pi{charge}_{x+1}_{q+1}_{z+1};1" not in key:
					continue
				if f"fMx_pi{charge}_{x+1}_{q+1}_{z+1};1" not in key:
					continue
				makePlot = True
				hist = inFile[f"hMx_pi{charge}_{x+1}_{q+1}_{z+1};1"]
				fFit = inFile_R.Get(f"fMx_pi{charge}_{x+1}_{q+1}_{z+1};1")
				
				fit = [fFit.GetParameter(i) for i in range(6)]
				values = hist.values()
				values[values == 0] = np.nan
				errors = hist.errors()
				errors[errors==0] = np.nan

				xEdges = hist.axis().edges()
				binCenters = (xEdges[:-1]+xEdges[1:])/2

				xVals = np.linspace(0, 1.45, 1000)
				#f_vals = [eval_func(fit, xVal) for xVal in xVals]

				axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
				axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

				axs[reg, ch].tick_params(which='both', width=2)
				axs[reg, ch].tick_params(which='major', length=7)
				axs[reg, ch].tick_params(which='minor', length=4)
				
				#not showing errorbars for visual clarity
				axs[reg, ch].errorbar(binCenters, values, errors, c=color ,drawstyle='steps-mid')
				#axs[reg, ch].plot(xVals, [eval_func_rho(fit, xVal) for xVal in xVals], c='r')
				#axs[reg, ch].plot(xVals, [eval_func_back(fit, xVal) for xVal in xVals], c='g')
			#plt.show()
		if makePlot:
			fig.savefig(outBase+f"/all_fits_Mx_xB_{x}_q2_{q}.pdf")
		plt.close(fig)
