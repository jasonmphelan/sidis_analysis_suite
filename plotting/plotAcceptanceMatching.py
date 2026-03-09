import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys

def eval_func(params, xVal):
	return params[0] + params[1]/xVal
def poly3(params, xVal):
	out = 0
	for i in range(4):
		out += params[i]*xVal**i
	return out

inName = sys.argv[1]
#inName = '../data/acceptance_matching/matchCut2D.root'
out_dir =sys.argv[2]
inFile = uproot.open(inName)
inFile_R = ROOT.TFile.Open(inName, "READ")
inFile_max = ROOT.TFile.Open("../data/acceptance_map/acceptanceMap_allE_final.root", "READ")

keys = inFile.keys()

nSec = 6
plot_max = True
#for key in keys: 
#	num_check = any(ch.isdigit() for ch in key)
if "Pi2K" in inName or "K2Pi" in inName:
	nSec = 1
	plot_max = False

print(nSec)

ch_string = ['pip', 'pim']
tit_string = [r"d(e, e'$\pi^{+}$)", "d(e, e'$\pi^{-}$)"]



for sec in range(nSec):
	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(24,8), layout='constrained')
	
	hist = np.empty( 2, dtype='object')
	fMax = np.empty( 2, dtype='object')
	fMin = np.empty( 2, dtype='object') 
	fAccMax = np.empty( 2, dtype='object') 


	max_params = np.empty( 2, dtype='object')
	min_params = np.empty( 2, dtype='object')
	acc_params = np.empty( 2, dtype='object')

	values = np.empty( 2, dtype='object')
	max_vals = np.empty( 2, dtype='object')
	min_vals = np.empty( 2, dtype='object')
	acc_max = np.empty( 2, dtype='object')

	xVals = np.linspace(1.25, 4.75, 1000)

	for i in range(2):
		if nSec == 1:
			print("max_0"+ch_string[i]+";1" )
			hist[i] = inFile[ f"hTheta_P_"+ch_string[i]+";1" ]
			fMax[i] = inFile_R.Get( "max_0_"+ch_string[i]+";1" )
			fMin[i] = inFile_R.Get( "min_0_"+ch_string[i]+";1" )
			fAccMax[i] = inFile_max.Get( f"f_{ch_string[i]}_{sec+1}_param_max")
		else:
			hist[i] = inFile[ f"hTheta_P_sec_{sec}_{ch_string[i]};1" ]
			fMax[i] = inFile_R.Get( f"max_{sec}_{ch_string[i]};1" )
			fMin[i] = inFile_R.Get( f"min_{sec}_{ch_string[i]};1" )
			fAccMax[i] = inFile_max.Get( f"f_{ch_string[i]}_{sec+1}_param_max")
		max_params[i] = [fMax[i].GetParameter(j) for j in range(2)]
		min_params[i] = [fMin[i].GetParameter(j) for j in range(2)]
		acc_params[i] = [fAccMax[i].GetParameter(j) for j in range(4)]

		xEdges = hist[i].axis(0).edges()
		yEdges = hist[i].axis(1).edges()
		values[i] = hist[i].values()
		values[i][values[i] == 0] = np.nan
		values[i] = np.reshape( values[i], (len(xEdges) - 1, len(yEdges) - 1) )
		max_vals[i] = [eval_func(max_params[i], xVal) for xVal in xVals]
		min_vals[i] = [eval_func(min_params[i], xVal) for xVal in xVals]
		acc_max[i] = [poly3(acc_params[i], xVal) for xVal in xVals]

	for i in range(2):


		axs[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
		axs[i].tick_params(which='both', width=2)
		axs[i].tick_params(which='major', length=7)
		axs[i].tick_params(which='minor', length=4)
		axs[i].tick_params(labelsize='18')
		
		axs[i].set_ylim( [5, 35] )

		axs[i].pcolormesh(xEdges, yEdges, values[i].T, cmap='winter', shading='auto')
		if( i == 1 ):
			axs[i].plot( xVals, max_vals[0], c='r', linestyle='-', linewidth=4, label = tit_string[0] + ' Cut')
			axs[i].plot( xVals, max_vals[1], c='orange', linestyle='-', linewidth=4, label = tit_string[1] + ' Cut')
		else:
			axs[i].plot( xVals, max_vals[0], c='r', linestyle='-', linewidth=4)
			axs[i].plot( xVals, max_vals[1], c='orange', linestyle='-', linewidth=4)
		if(nSec!=1):
			axs[i].plot( xVals, acc_max[0], c='r', linestyle='--', linewidth=4)
			axs[i].plot( xVals, acc_max[1], c='orange', linestyle='--', linewidth=4)
		axs[i].plot( xVals, min_vals[0], c='r', linestyle='-', linewidth=4) 
		axs[i].plot( xVals, min_vals[1], c='orange', linestyle='-', linewidth=4) 

		lower = np.minimum(max_vals[0], max_vals[1])
		lower = np.minimum(lower, acc_max[0])
		lower = np.minimum(lower, acc_max[1])
		
		upper = np.maximum(min_vals[0], min_vals[1])

		if(nSec != 1 ):
			axs[i].fill_between(
				xVals,
				lower,
				upper,
				color='magenta',
				alpha=0.3,
				linewidth=0
			)


		axs[i].set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=24)
		axs[i].set_ylabel(r"$\theta_{\pi}$ [deg.]", fontsize=24)
		if(nSec != 1):
			axs[i].set_title(f"Sector {sec}, " + tit_string[i], fontsize=22)
		else:
			axs[i].set_title(tit_string[i], fontsize=22)
		
	fig.legend(fontsize=24)

	if nSec == 1:
		if "Pi2K" in inName:
			fig.savefig(f'{out_dir}/hMatch_pions.png')
		else:
			fig.savefig(f'{out_dir}/hMatch_kaons.png')
	else:
		fig.savefig(f'{out_dir}/hMatch_{sec+1}.png')
