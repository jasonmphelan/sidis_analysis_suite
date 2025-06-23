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
inFile = uproot.open(sys.argv[2])

colorList = ['red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

for q in range(12):
	
	fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
	ax[1].set_ylim( [0.8, 1.2] )
	ax[0].set_ylim( [0, 1] )
	ax[1].set_xlim( [.3, .8] )
		
	z = np.linspace( .3, 1, 500 )
	ax[0].plot( z, np.asarray(ff(z)), color = 'black', linestyle = '--')
	ax[1].axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')
	
	ax[1].set_xlabel(r'$z$')
	ax[1].set_ylabel(rf'$r(z, x_B)/$FF($z$))')
	ax[0].set_ylabel(r'$r(z)$')

	ax[0].text(.7, .7, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*(q), 2+.5*(q+1)), fontsize=12)
	
	for q in range(12):
	

		hist = inFile[f'hRatio_{q+1}_{x+1}']
		
		values =np.array( hist.values() )
		values[values <= 0] = np.nan
		if( np.isnan(values).all() ):
			continue
		errors = np.array(hist.errors())

		xEdges = hist.axis().edges()
		binCenters = np.asarray((xEdges[:-1]+xEdges[1:])/2 + ( -.01 + .002*q))


		label_x = r'%.2f $< x_B <$ %.2f'%( .01 + .04*x, .1 + .04*(x+1))


		ax[0].errorbar( binCenters, values, errors, label=label_x, marker='o',
						color = colorList[x], linestyle = '',capsize = 2, lw = 1, capthick = 1)

		ratio = np.asarray(ff(binCenters))/values
		ratio_err =  ratio*np.sqrt( (errors/values)**2)
		ratio_err[ratio_err<=0] = np.nan

		ax[1].errorbar( binCenters, ratio, ratio_err, marker='o',
							color = colorList[x], linestyle = '',capsize = 2, lw = 1, capthick = 1)

			

		#plt.errorbar(binCenters, ratio, yerr=ratio_err, color=colorList[x], label=tit, linestyle='',marker='.', markersize=10, mec='black', capsize=2)
		#plt.plot(binCenters, ratio, color=colorList[q], label=tit, linestyle='',marker='.', markersize=10, mec='black')
		

	ax[0].legend()
		#plt.text(.95, .5, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1) ), fontsize=12)
		#plt.title( r'%.1f $< Q^2 <$ %.1f [GeV$^2$]'%( 2 + 0.5*q, 2 + 0.5*(q+1) ), fontsize=16)
	fig.savefig('{0}/ratio_{1}.pdf'.format(outFileName, q+1))
		
	plt.close()
