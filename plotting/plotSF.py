import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
inFile = uproot.open("../data/SF_fits.root")
inFile_R = ROOT.TFile.Open("../data/SF_fits.root", "READ")

obj_base = "hSF_"
def eval_func(params, xVal):
	return params[0] + params[1]*xVal + params[2]*xVal**2

fig, axs = plt.subplots(3, 2, figsize=(14,16), layout='constrained')
for sec in range(6):
	reg = sec%3
	ch = math.floor(sec/3)
	
	axs[reg,ch].set_ylim( [0.1, .35] )
	#plt.subplots_adjust( wspace = 0.25, hspace=.25)
	hSF = inFile[f"hSF_sec_{sec}"]
	hMean = inFile_R.Get(f"fMean_SF_{sec}")
	hSig = inFile_R.Get(f"fSigma_SF_{sec}")
	
	mean_params = [hMean.GetParameter(i) for i in range(3)]
	sig_params = [hSig.GetParameter(i) for i in range(3)]

	values = hSF.values()
	values[values == 0] = np.nan

	xEdges = hSF.axis(0).edges()
	yEdges = hSF.axis(1).edges()

	xVals = np.linspace(.8, xEdges[-1], 1000)
	mean_vals = [eval_func(mean_params, xVal) for xVal in xVals]
	max_vals = [3.5*eval_func(sig_params, xVal) + eval_func(mean_params, xVal) for xVal in xVals]
	min_vals = [-3.5*eval_func(sig_params, xVal) + eval_func(mean_params, xVal) for xVal in xVals]

	values = np.reshape( values, (len(xEdges) - 1, len(yEdges) - 1) )

	axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
	axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())

	axs[reg, ch].tick_params(which='both', width=2)
	axs[reg, ch].tick_params(which='major', length=7)
	axs[reg, ch].tick_params(which='minor', length=4)
	
	axs[reg, ch].pcolormesh(xEdges, yEdges, values.T, shading='auto')
	axs[reg, ch].plot( xVals, mean_vals, c='r', linestyle='-', linewidth=2)
	axs[reg, ch].plot( xVals, max_vals, c='r', linestyle='-', linewidth=2)
	axs[reg, ch].plot( xVals, min_vals, c='r', linestyle='-', linewidth=2)
	
	axs[reg, ch].set_ylabel(r"SF", fontsize=14)
	axs[reg, ch].set_xlabel(r"$p_{e}$ [GeV]", fontsize=14)
	#if ch != 0:
	#	axs[reg, ch].set_ylabel("")
	#	axs[reg, ch].set_yticks(color='w')
	axs[reg, ch].set_title(f"Sector {sec + 1}")

#plt.show()
fig.savefig(f"selection_plots/hSF.pdf")
