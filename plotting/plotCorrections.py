import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as LogNorm
import matplotlib.ticker as ticker
import math
import sys 
#tools
def getTitle(var):

	titdic = {
		'Q2':r'$Q^2$ [GeV$^2$]',
		'xB':r'$x_{B}$',
		'Z':'$z$',
		'M_x':'$M_{X}$ [GeV$^2$]',
		'omega':r'$\omega$',
		'W2':r'$W^2$',
		'y':r'$y$',
		'p_e':r'$|p_{e}|$ [GeV]',
		'p_pi':r'$|p_{\pi}|$ [GeV]',
		'w-':r'$w^-$',
		'w+':r'$w^+$',
		'w+/w-':r'$w^+/w^-$',
		}
	return titdic[var]


def ff(z):
	return (1-z)/(z+1)

def makePlots( corrType, hist, ext):
	

	values = hist.values()
	errors = hist.errors()

	values[values == 0] = np.nan
	errors[errors == 0] = np.nan
	
	xBins = 14 
	if 'Kaon' in hist.name:
		xBins=7
	qBins = 12 

	zEdges = hist.axis(2).edges()
	zBinCenters = (zEdges[:-1]+zEdges[1:])/2

	textHeight = 0.4

	outDir = ''
	textHeight = 0
	yMin = 0
	yMax = 10

	deltaX = 0.04

	wType = ''
	if hist.name[-1] == 'P' or hist.name[-3] == 'P':
		wType = 'w+'
	elif hist.name[-1] == 'M' or hist.name[-3] == 'M':
		wType = 'w-'
	elif 'Kaon' in hist.name:
		wType = 'w+/w-'
		yMin = 0
		yMax = 1.1	
	else:
		wType = 'w+/w-'
		yMin = .5
		yMax = 1.5

	if corrType == 'bin':
		outDir = 'bin_migration_corrections/'
		textHeight = .75
		yMin = .5
		yMax = 1.5
	if corrType == 'acc' or corrType == 'mc':
		outDir = 'acceptance_corrections/'
		textHeight = 5 
		if wType == 'w+/w-':
			textHeight = 1.2
	if corrType == 'k2pi':
		outDir = 'k_to_pi_corrections/'
		textHeight = .75
		yMin = 0
		yMax = 1.05
		deltaX = 0.08
	if corrType == 'pi2k':
		outDir = 'pi_to_k_corrections/'
		textHeight = .75
		yMin = 0.5
		yMax = 1.05
		deltaX = 0.08
	if 'ratio' in ext:
		yMin = 0.75
		yMax = 1.25
		textHeight = 1.2
	outDir = f'/volatile/clas12/users/jphelan/SIDIS/analysis_note/corrections_{corrType}/'
    

	fig, axs = plt.subplots( 4, 3,figsize=(18, 10))
	plt.subplots_adjust( wspace=.1, hspace=.1 )
	
	for q in range(qBins):
		for x in range(xBins):
			axs[math.floor(q/3), q%3].text(.75, textHeight, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*(q), 2+.5*(q+1)), fontsize=12)
				
			axs[math.floor(q/3), q%3].plot( zBinCenters, values[x][q], colorList[x], marker = '.', linestyle='none', label = r'%0.2f $< x_B <$ %0.2f'%(.1+deltaX*x, .1 + deltaX*(x+1)), ms=10, mec='black' )
			axs[math.floor(q/3), q%3].errorbar(zBinCenters, values[x][q], yerr=errors[x][q], color = colorList[x], linestyle = '',capsize = 2, lw = 1, capthick = 1)

			if corrType == 'k2pi' or corrType == 'pi2k':
				bin_num = int(hist.name[-1])
				axs[math.floor(q/3), q%3].text(.745, textHeight - .1, r'%.2f $< p_{\pi} <$ %.2f'%( p_bin[bin_num], p_bin[bin_num+1]), fontsize=12)
			
			axs[math.floor(q/3), q%3].axhline( y = 1, color = 'black', linestyle = '--')
	

	

	for ax in axs.flat:
		ax.set_xlabel( getTitle('Z'), fontsize=16 )
		ax.set_ylabel( getTitle(wType), fontsize=16)
		ax.label_outer()

		ax.set_ylim( [yMin, yMax] )
		ax.set_xlim( [0.3, 1] )

	axs[0,1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.5), ncol=7, fancybox=True, shadow=True)	
	print('Writing : ' + outDir+hist.name+ext)	
	plt.savefig(outDir+hist.name+ext, bbox_inches='tight')
	
		

#colorList = ['red', 'blue', 'magenta', 'green', 'gold', 'cyan', 'black', 'blueviolet','darkorange', 'brown','red', 'blue', 'magenta','green']
colorList = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']

p_bin = [1.25, 2.00, 3, 4, 5.00]
##############################################
#start script#
#############################################

inFile_name = sys.argv[1]#input('File name: \n')
corrType = sys.argv[2]#input('Correction type (bin, acc, k2pi, pi2k): \n')

inFile = uproot.open('../data/correctionFiles/'+inFile_name)
print('../data/correctionFiles/'+inFile_name)
keyList = []

if corrType == 'bin':
	keyList = ['hBinMigration', 'hBinMigrationP', 'hBinMigrationM']
if corrType == 'acc':
	keyList = ['hAccCorrection', 'hAccCorrectionP', 'hAccCorrectionM']
if corrType == 'k2pi' or corrType == 'pi2k':
	keyList = ['hKaonCorr', 'hKaonCorrP', 'hKaonCorrM']
if corrType == 'mc':
	keyList = ['hMcCorrection', 'hMcCorrectionP', 'hMcCorrectionM']

energy = '_'
if '3d' in inFile_name:
	energy = energy + '3d_'
if 'no_match' in inFile_name:
	energy = energy + 'no_match.pdf'
if 'ratio' in inFile_name:
	energy = energy + 'ratio.pdf'
if '10.2' in inFile_name:
	energy = energy + '10.2.pdf'
elif '10.4' in inFile_name:
	energy = energy + '10.4.pdf'
elif '10.6' in inFile_name:
	energy = energy + '10.6.pdf'
if energy == '_':
	energy = '.pdf'

print( energy )

for key in keyList:
	if 'k2pi' in inFile_name or 'pi2k' in inFile_name:
		for p in range(4):
			hist = inFile[key+f'_{p};1']
			makePlots( corrType, hist, energy)
			
	else:
		hist = inFile[key]
		makePlots( corrType, hist, energy)
