import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import csv 

def ff(z):
	return (1-z)/(z+1)

inFile = uproot.open("../sidis_plotter/rootFiles/ratios/ratios_10.2_data_kaon_corrected.root")
inFile_rho = uproot.open("../sidis_plotter/rootFiles/ratios/ratios_10.2_data_rho_corrected.root")

outFileName = 'rho_effect_xb'
colorList = ['red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

for q in range(12):
	plt.figure(figsize=(12,4))
	ax = plt.gca()
	ax.set_ylim( [0.8, 1.2] )
	ax.set_xlim( [.3, .8] )
		
	z = np.linspace( .3, 1., 500 )
	plt.axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')
	plt.xlabel(r'$z$')
	plt.ylabel(r'$r$(MC + Kaon + Rho)/$r$(MC + Kaon)')
	
	for x in range(14):
	
		hist = inFile[f'hRatio_{q+1}_{x+1}']
		hist_rho = inFile_rho[f'hRatio_{q+1}_{x+1}']
		
		values =np.array( hist.values() )
		values[values <= 0] = np.nan
		errors = np.array(hist.errors())

		values_rho = np.array(hist_rho.values())
		values_rho[values_rho <= 0] = np.nan
		errors_rho = np.array(hist_rho.errors())
		
		
		ratio = values_rho/values
		ratio_err =  ratio*np.sqrt( (errors/values)**2 + (errors_rho/values_rho)**2 )
		ratio_err[ratio_err<=0] = np.nan
	
		xEdges = hist.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2 + ( -.01 + .002*q)
	
		if( np.isnan(ratio).all() ):
			continue
		tit = r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1) )

		plt.errorbar(binCenters, ratio, yerr=ratio_err, color=colorList[x], label=tit, linestyle='',marker='.', markersize=10, mec='black', capsize=2)
		#plt.plot(binCenters, ratio, color=colorList[q], label=tit, linestyle='',marker='.', markersize=10, mec='black')
		

	plt.legend()
	#plt.text(.95, .5, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1) ), fontsize=12)
	plt.title( r'%.1f $< Q^2 <$ %.1f [GeV$^2$]'%( 2 + 0.5*q, 2 + 0.5*(q+1) ), fontsize=16)
	plt.savefig('{0}/ratio_{1}.pdf'.format(outFileName, q+1))
		
	plt.close()
