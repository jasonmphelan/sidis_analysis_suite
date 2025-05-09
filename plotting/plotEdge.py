import uproot
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import sys 
in_name = int(sys.argv[1])
inFile = uproot.open(in_name)
out_dir = int(sys.argv[2])
#inFile = uproot.open("../histograms/analysis_note/fiducials.root.root")

fid_cuts = [[5, 5, 10], [4, 4, 10], [4, 4, 10] ]

parList = ["e", "pi", "pi"]
stList = ["_pip", "_pip", "_pim"]
labList=['All', r'$5<\theta <12$',r'$12<\theta <18$', r'$18<\theta < 24$', r'$24<\theta<40$']
colorList = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']

for par in range(3):
	#plt.subplots_adjust( wspace = 0.25, hspace=.25)
	for reg in range(3):
		fig, axs = plt.subplots(2, 3, figsize=(12,10), layout='constrained')
		for sec in range(6):
			yPl = sec%3
			xPl = math.floor(sec/3)
			for tBin in range(5):

				hEdge = inFile[f"hFid_{parList[par]}_w_reg_{reg}_{tBin}_{sec}{stList[par]}"]

				values = hEdge.values()
				errors = hEdge.errors()

				xEdges = hEdge.axis().edges()
				binCenters = (xEdges[:-1]+xEdges[1:])/2
		

				#plt.figure(figsize=(8,6))
				axs[xPl, yPl].xaxis.set_minor_locator(ticker.AutoMinorLocator())
				axs[xPl, yPl].yaxis.set_minor_locator(ticker.AutoMinorLocator())
		
				axs[xPl, yPl].tick_params(which='both', width=2)
				axs[xPl, yPl].tick_params(which='major', length=7)
				axs[xPl, yPl].tick_params(which='minor', length=4)
				
				axs[xPl, yPl].errorbar(binCenters, values, errors, capsize=1, fmt=colorList[tBin], ecolor=colorList[tBin], linestyle='', marker='o',markersize=1, label=labList[tBin])
				
				axs[xPl, yPl].set_xlabel(r"Edge [mm]")
				axs[xPl, yPl].set_ylabel(r"$\langle \chi^2 /NDF \rangle$")
				#if ch != 0:
				#	axs[reg, ch].set_ylabel("")
				#	axs[reg, ch].set_yticks(color='w')
				axs[xPl, yPl].set_title(rf"Sector {sec + 1}, Region {reg+1}")
		
			axs[xPl, yPl].axvline( x=fid_cuts[par][reg], ymin=0, ymax=1, color='k', linewidth=2, linestyle='-' )

		plt.margins(y=0)
		plt.margins(x=0)

		#plt.show()
		plt.legend()
		fig.savefig(f"{out_dir}/hFid_{par}{stList[par]}_reg_{reg}.pdf")
