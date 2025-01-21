import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math


def eval_func(params, xVal):
	return params[0] + params[1]/xVal
inName = '../data/acceptance_matching/matchCut2D.root'
inFile = uproot.open(inName)
inFile_R = ROOT.TFile.Open(inName, "READ")

keys = inFile.keys()

nSec = 6

for key in keys: 
	num_check = any(ch.isdigit() for ch in key)
	if "Kaon" in inName or "kaon" in inName:
		nSec = 1
		break


ch_string = ['pip', 'pim']
tit_string = [r"d(e, e'$\pi^{+}$)", "d(e, e'$\pi^{-}$)"]



for sec in range(nSec):
	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(24,8), layout='constrained')
	
	hist = np.empty( 2, dtype='object')
	fMax = np.empty( 2, dtype='object')
	fMin = np.empty( 2, dtype='object') 

	max_params = np.empty( 2, dtype='object')
	min_params = np.empty( 2, dtype='object')
	
	values = np.empty( 2, dtype='object')
	max_vals = np.empty( 2, dtype='object')
	min_vals = np.empty( 2, dtype='object')
	xVals = np.linspace(1.25, 4.75, 1000)

	for i in range(2):
		if nSec == 1:
			hist[i] = inFile[ f"hTheta_P_"+ch_string[i] ]
			fMax[i] = inFile_R.Get( "max"+ch_string[i] )
			fMin[i] = inFile_R.Get( "min"+ch_string[i] )
		else:
			hist[i] = inFile[ f"hTheta_P_sec_{sec}_{ch_string[i]}" ]
			fMax[i] = inFile_R.Get( f"max_{sec}_{ch_string[i]}" )
			fMin[i] = inFile_R.Get( f"min_{sec}_{ch_string[i]}" )
		max_params[i] = [fMax[i].GetParameter(j) for j in range(2)]
		min_params[i] = [fMin[i].GetParameter(j) for j in range(2)]
		xEdges = hist[i].axis(0).edges()
		yEdges = hist[i].axis(1).edges()
		values[i] = hist[i].values()
		values[i][values[i] == 0] = np.nan
		values[i] = np.reshape( values[i], (len(xEdges) - 1, len(yEdges) - 1) )
		max_vals[i] = [eval_func(max_params[i], xVal) for xVal in xVals]
		min_vals[i] = [eval_func(min_params[i], xVal) for xVal in xVals]

	for i in range(2):


		axs[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[i].tick_params(which='both', width=2)
		axs[i].tick_params(which='major', length=7)
		axs[i].tick_params(which='minor', length=4)

		axs[i].set_ylim( [5, 35] )
		axs[i].pcolormesh(xEdges, yEdges, values[i].T, cmap='winter', shading='auto')
		if( i == 1 ):
			axs[i].plot( xVals, max_vals[0], c='r', linestyle='-', linewidth=4, label = tit_string[0] + ' Cut')
			axs[i].plot( xVals, max_vals[1], c='orange', linestyle='-', linewidth=4, label = tit_string[1] + ' Cut')
		else:
			axs[i].plot( xVals, max_vals[0], c='r', linestyle='-', linewidth=4)
			axs[i].plot( xVals, max_vals[1], c='orange', linestyle='-', linewidth=4)
		axs[i].plot( xVals, min_vals[0], c='r', linestyle='-', linewidth=4) 
		axs[i].plot( xVals, min_vals[1], c='orange', linestyle='-', linewidth=4) 


		axs[i].set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=16)
		axs[i].set_ylabel(r"$\theta_{\pi}$ [deg.]", fontsize=16)
		if(nSec != 1):
			axs[i].set_title(f"Sector {sec}, " + tit_string[i], fontsize=16)
		else:
			axs[i].set_title(tit_string[i], fontsize=16)
		
	fig.legend(fontsize=16)

	if nSec == 1:
		fig.savefig(f'acceptance_matching/hMatch_kaons.pdf')
	else:
		fig.savefig(f'acceptance_matching/hMatch_{sec+1}.pdf')
