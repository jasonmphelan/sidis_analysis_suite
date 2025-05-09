import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import sys 

in_name = sys.argv[1]
out_dir = sys.argv[2]
#inFile = uproot.open("../histograms/analysis_note/detector_plots_10.4.root")
inFile = uproot.open(in_name)

parList = ["e", "pi"]
chList = ["pip", "pim"]
stList = ["\pi^+", "\pi^-"]

for par in parList:
	fig, axs = plt.subplots(3, 2, figsize=(10,12), layout='constrained')
	#plt.subplots_adjust( wspace = 0.25, hspace=.25)
	for reg in range(3):
		for ch in range(2):
			hFid_bef_reg_0 = inFile[f"hFid_{par}_bef_reg_{reg}_{chList[ch]}"]
			hFid_aft_reg_0 = inFile[f"hFid_{par}_aft_reg_{reg}_{chList[ch]}"]

			values_bef = hFid_bef_reg_0.values()
			values_bef[values_bef == 0] = np.nan
			values_aft = hFid_aft_reg_0.values()
			values_aft[values_aft == 0] = np.nan

			xEdges = hFid_bef_reg_0.axis(0).edges()
			yEdges = hFid_bef_reg_0.axis(1).edges()

			values_bef = np.reshape( values_bef, (len(xEdges) - 1, len(yEdges) - 1) )

			#plt.figure(figsize=(8,6))
			axs[reg, ch].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			axs[reg, ch].yaxis.set_minor_locator(ticker.AutoMinorLocator())
	
			axs[reg, ch].tick_params(which='both', width=2)
			axs[reg, ch].tick_params(which='major', length=7)
			axs[reg, ch].tick_params(which='minor', length=4)
			
			axs[reg, ch].pcolormesh(xEdges, yEdges, values_bef.T, cmap='autumn', shading='auto')
			axs[reg, ch].pcolormesh(xEdges, yEdges, values_aft.T, cmap='viridis', shading='auto', norm="log")
			axs[reg, ch].set_xlabel(r"X$_{DC}$ [mm]")
			axs[reg, ch].set_ylabel(r"Y$_{DC}$ [mm]")
			#if ch != 0:
			#	axs[reg, ch].set_ylabel("")
			#	axs[reg, ch].set_yticks(color='w')
			axs[reg, ch].set_title(rf"$(e, e'{stList[ch]})$, Region {reg+1}")

	#plt.show()
	fig.savefig(f"{out_dir}/hFid_{par}.pdf")
