import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import sys 

in_name = (sys.argv[1])
out_dir = (sys.argv[2])

def setTitle2D( ax, key):
	if "Beta" in key or "Min" in key or "Pos" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]", fontsize=18)
		ax.set_ylabel(r"$\beta$", fontsize=18)
def setTitle1D( ax, key):
	ax.set_ylabel("Counts [a.u.]", fontsize=18)
	if "Beta" in key or "Min" in key or "Pos" in key:
		ax.set_xlabel(r"$\beta$", fontsize=18)
		ax.set_yscale("log")

inFile = uproot.open(in_name)
for key in inFile.keys():
	if "pim" in key:
		continue
	if "_e_" not in key:
		continue
	if ("Beta" not in key):
		continue
	if ";1" in key:# or "Min" not in key or "Pos" not in key:
		continue
	print(key)
	hist = inFile[key]
	if isinstance( hist, uproot.models.TH.Model_TH2F_v4): 
		if "reg" in key:
			continue
		
		fig, axs = plt.subplots(2, 3, figsize=(18,12), layout='constrained')
		
		eName = key
		pipName = eName.replace("_e_","_pip_")
		pimName = eName.replace("_e_","_pim_")

	
		keys = [eName, pipName, pimName]
		labels = ['e', '\pi^+', '\pi^-']

		for i in range(3):
			hist2d = inFile[keys[i]]
			values = hist2d.values()
			values[values == 0] = np.nan

			name_1 = keys[i].replace("_p_", "_")
			hist1d = inFile[name_1.replace(';2', ';1')]

			xEdges = hist2d.axis(0).edges()
			yEdges = hist2d.axis(1).edges()
		
			axs[0, i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
			axs[0,i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
			
	
			axs[0,i].tick_params(which='both', width=2)
			axs[0,i].tick_params(which='major', length=7)
			axs[0,i].tick_params(which='minor', length=4)
		
		
			axs[0,i].pcolormesh(xEdges, yEdges, values.T, cmap='viridis', shading='auto', norm="log")
			#axs[i].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
			axs[0,i].set_xlabel(fr"$p_{labels[i]}$ [GeV]", fontsize=18)
			axs[0,i].set_ylabel(fr"$\beta_{labels[i]}$", fontsize=18)


			values_1 = hist1d.values()
			errors_1 = hist1d.errors()
		

			xEdges_1 = hist1d.axis().edges()
			binCenters = (xEdges_1[:-1]+xEdges_1[1:])/2
		
			axs[1,i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
	
			axs[1,i].tick_params(which='both', width=2, labelsize=12)
			axs[1,i].tick_params(which='major', length=7, labelsize=12)
			axs[1,i].tick_params(which='minor', length=4, labelsize=12)
		

			#axs[0].plot(binCenters, values_pip, 'k', drawstyle='steps-mid')
			axs[1,i].errorbar(binCenters, values_1, errors_1, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid')
			axs[1,i].set_ylabel("Counts", fontsize=18)
			axs[1,i].set_xlabel(fr"$\beta_{labels[i]}$", fontsize=18)
			axs[1,i].set_yscale("log")

		pdfName = eName.replace(";2", "")
		pdfName = pdfName.replace("_pip", "")
		print(f"{out_dir}/{pdfName}.pdf")
		fig.savefig(f"{out_dir}/{pdfName}.pdf")
		plt.close(fig)
		
'''
fig, axs = plt.subplots(2, 3, figsize=(18,12), layout='constrained')

eName = 'hEl_b_p;2'
pipName = 'hPos_b_p;2'
pimName = 'hMin_b_p;2'


keys = [eName, pipName, pimName]
labels = ['_e', '^+', '^-']

for i in range(3):
	hist2d = inFile[keys[i]]
	values = hist2d.values()
	values[values == 0] = np.nan

	name_1 = keys[i].replace("_p_", "_")
	hist1d = inFile[name_1.replace(';2', ';1')]

	xEdges = hist2d.axis(0).edges()
	yEdges = hist2d.axis(1).edges()

	axs[0, i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
	axs[0,i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
	

	axs[0,i].tick_params(which='both', width=2)
	axs[0,i].tick_params(which='major', length=7)
	axs[0,i].tick_params(which='minor', length=4)


	axs[0,i].pcolormesh(xEdges, yEdges, values.T, cmap='viridis', shading='auto', norm="log")
	#axs[i].set_title(rf"$d(e, e'\pi^+)$", fontsize=16)
	axs[0,i].set_xlabel(fr"$p{labels[i]}$ [GeV]", fontsize=18)
	axs[0,i].set_ylabel(fr"$\beta{labels[i]}$", fontsize=18)


	values_1 = hist1d.values()
	errors_1 = hist1d.errors()


	xEdges_1 = hist1d.axis().edges()
	binCenters = (xEdges_1[:-1]+xEdges_1[1:])/2

	axs[1,i].xaxis.set_minor_locator(ticker.AutoMinorLocator())

	axs[1,i].tick_params(which='both', width=2, labelsize=12)
	axs[1,i].tick_params(which='major', length=7, labelsize=12)
	axs[1,i].tick_params(which='minor', length=4, labelsize=12)


	#axs[0].plot(binCenters, values_pip, 'k', drawstyle='steps-mid')
	axs[1,i].errorbar(binCenters, values_1, errors_1, capsize=2, fmt='k', ecolor='k',drawstyle='steps-mid')
	axs[1,i].set_ylabel("Counts", fontsize=18)
	axs[1,i].set_xlabel(fr"$\beta{labels[i]}$", fontsize=18)
	axs[1,i].set_yscale("log")

pdfName = 'hBeta_p_All_tracks'
print(f"{out_dir}/{pdfName}.pdf")
fig.savefig(f"{out_dir}/{pdfName}.pdf")
plt.close(fig)

'''

