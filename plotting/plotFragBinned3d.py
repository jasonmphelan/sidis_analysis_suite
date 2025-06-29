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

def getVar(name):
	if name in ['pT', 'Pt', 'pt', 'PT']:
		return r'$p_T$', 4, 0, 1.2, 'GeV'
	if name in ['pi_sector', 'sector_pi']:
		return '', 6, 0, 6, ''
	if name in ['e_sector', 'sector_e']:
		return '', 6, 0, 6, ''

#python makeRatios.py outfile file_1 title_1 file_2 title_2 ...
outFileName = sys.argv[1]
inFile = uproot.open(sys.argv[2])
var = sys.argv[3]

lab, nBins, min, max, unit = getVar(var)
delta = (max - min)/nBins
colorList = ['red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

for x in range(14):
	for q in range(12):
		fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		ax[1].set_ylim( [0.8, 1.2] )
		ax[0].set_ylim( [0, 1] )
		ax[1].set_xlim( [.3, .8] )
		
		z = np.linspace( .3, 1, 500 )
		ax[0].plot( z, np.asarray(ff(z)), color = 'black', linestyle = '--')
		ax[1].axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')
	
		ax[1].set_xlabel(r'$z$')
		ax[1].set_ylabel(rf'$r(z, Q^2)/$FF(z))')
		ax[0].set_ylabel(r'$r(z)$')
	

	
		ax[0].text(.7, .7, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1)), fontsize=12)
		ax[0].text(.7, .65, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*(q), 2+.5*(q+1)), fontsize=12)
		for bin in range(nBins):

			hist = inFile[f'hRatio_{q+1}_{x+1}_{bin+1}']
		
			values =np.array( hist.values() )
			values[values <= 0] = np.nan
			if( np.isnan(values).all() ):
				continue
			errors = np.array(hist.errors())

			xEdges = hist.axis().edges()
			binCenters = np.asarray((xEdges[:-1]+xEdges[1:])/2)

			label = ''

			
			if 'sector' in var:
				if 'pi' in var:
					label = rf'$\pi$ Sector {bin+1}'
				elif 'e' in var:
					label = rf'$e$ Sector {bin+1}'
				else:
					label = rf'Sector {bin+1}'
			else:
				label = rf'${bin + delta*bin} < ${lab} $< {bin + delta*(bin+1)}$ {unit}'

			ax[0].errorbar( binCenters, values, errors, marker='o', label=label,
						color = colorList[q], linestyle = '',capsize = 2, lw = 1, capthick = 1)

			ratio = np.asarray(ff(binCenters))/values
			ratio_err =  ratio*np.sqrt( (errors/values)**2)
			ratio_err[ratio_err<=0] = np.nan

			ax[1].errorbar( binCenters, ratio, ratio_err, marker='o',
							color = colorList[q], linestyle = '',capsize = 2, lw = 1, capthick = 1)

		

    	ax[0].legend()
		fig.savefig('{0}/ratio_{1}_{2}.pdf'.format(outFileName, q+1, x+1))
		
		plt.close()
