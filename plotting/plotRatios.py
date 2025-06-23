import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import csv 
import sys


def ff(z):
	return (1-z)/(z+1)

#python makeRatios.py outfile file_1 title_1 file_2 title_2 ...
outFileName = sys.argv[1]
nFiles = int((len(sys.argv) - 2)/2)

inFile_list = [sys.argv[2*i + 2] for i in range(nFiles)]
inFile_tits = [sys.argv[2*i + 3] for i in range(nFiles)]


colorList = ['red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

for q in range(12):
	
	for x in range(14):
	
		fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		ax[1].set_ylim( [0.8, 1.2] )
		ax[0].set_ylim( [0, 1] )
		ax[1].set_xlim( [.3, .8] )
		
		z = np.linspace( .3, 1, 500 )
		ax[0].plot( z, np.asarray(ff(z)), color = 'black', linestyle = '--')
		ax[1].axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')
	
		ax[1].set_xlabel(r'$z$')
		ax[1].set_ylabel(rf'$r(z)/r$({inFile_tits[0]})')
		ax[0].set_ylabel(r'$r(z)$')

		ax[0].text(.7, .75, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*(q), 2+.5*(q+1)), fontsize=12)
		ax[0].text(.7, .7, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1)), fontsize=12)

		inFile = uproot.open(inFile_list[0])
		hist = inFile[f'hRatio_{q+1}_{x+1}']
		
		values =np.array( hist.values() )
		values[values <= 0] = np.nan
		errors = np.array(hist.errors())

		xEdges = hist.axis().edges()
		binCenters = (xEdges[:-1]+xEdges[1:])/2 + ( -.01 + .002*q)

		ax[0].errorbar( binCenters, values, errors, label=inFile_tits[0], marker='o',
						color = colorList[0], linestyle = '',capsize = 2, lw = 1, capthick = 1)

		if( np.isnan(values).all() ):
			continue

		for num in range( 1, len(inFile_list) ):
			inFile_temp = uproot.open(inFile_list[num])
			hist_temp = inFile_temp[f'hRatio_{q+1}_{x+1}']
		
			values_temp =np.array( hist_temp.values() )
			values_temp[values_temp <= 0] = np.nan
			errors_temp = np.array(hist_temp.errors())
			ratio = values_temp/values
			ratio_err =  ratio*np.sqrt( (errors/values)**2 + (errors_temp/values_temp)**2 )
			ratio_err[ratio_err<=0] = np.nan
				
			ax[0].errorbar( binCenters, values_temp, errors_temp, label=inFile_tits[num], marker='o',
							color = colorList[num], linestyle = '',capsize = 2, lw = 1, capthick = 1)

			ax[1].errorbar( binCenters, ratio, ratio_err, marker='o',
							color = colorList[num], linestyle = '',capsize = 2, lw = 1, capthick = 1)

			

		#plt.errorbar(binCenters, ratio, yerr=ratio_err, color=colorList[x], label=tit, linestyle='',marker='.', markersize=10, mec='black', capsize=2)
		#plt.plot(binCenters, ratio, color=colorList[q], label=tit, linestyle='',marker='.', markersize=10, mec='black')
		

		ax[0].legend()
		#plt.text(.95, .5, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1) ), fontsize=12)
		#plt.title( r'%.1f $< Q^2 <$ %.1f [GeV$^2$]'%( 2 + 0.5*q, 2 + 0.5*(q+1) ), fontsize=16)
		fig.savefig('{0}/ratio_{1}_{2}.pdf'.format(outFileName, q+1, x+1))
		
		plt.close()
