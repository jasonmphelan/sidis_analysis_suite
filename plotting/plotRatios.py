import uproot
import numpy as np
import matplotlib.pyplot as plt
import sys


def ff(z):
	return (1-z)/(z+1)

def ratio_to_ff(R, R_err):
	"""Convert multiplicity ratio R = pip/pim to fragmentation function (4-R)/(4R-1)."""
	FF = (4 - R) / (4*R - 1)
	# dFF/dR = -15 / (4R-1)^2
	denom = 4*R - 1
	FF_err = np.abs(15 / denom**2) * R_err
	return FF, FF_err

def read_hist(inFile, key):
	"""Return (values, errors) arrays for a histogram key, or (None, None) if absent."""
	if key not in inFile.keys():
		return None, None
	h = inFile[key]
	return np.array(h.values()), np.array(h.errors())

def read_xs(inFile, q, x):
	"""Read pip and pim cross section histograms, apply kaon and rho corrections,
	and return FF with errors. Kaon histograms are added; rho histograms are subtracted.
	Both corrections are applied only when the histograms are present."""
	key_pip = f'hRatio_{q+1}_{x+1};1'
	key_pim = f'hRatioPim_{q+1}_{x+1};1'

	pip, pip_err = read_hist(inFile, key_pip)
	pim, pim_err = read_hist(inFile, key_pim)
	if pip is None or pim is None:
		return None

	xEdges = inFile[key_pip].axis().edges()
	binCenters = (xEdges[:-1] + xEdges[1:]) / 2

	# Add kaon corrections (independent uncertainties)
	k_pip, k_pip_err = read_hist(inFile, f'hRatio_k_{q+1}_{x+1};1')
	k_pim, k_pim_err = read_hist(inFile, f'hRatioPim_k_{q+1}_{x+1};1')
	if k_pip is not None:
		pip     = pip     + k_pip
		pip_err = np.sqrt(pip_err**2 + k_pip_err**2)
	if k_pim is not None:
		pim     = pim     + k_pim
		pim_err = np.sqrt(pim_err**2 + k_pim_err**2)

	# Subtract rho background — same histogram for pip and pim, so correlated.
	# R = (pip - r) / (pim - r)
	# dR/dr = (R - 1) / pim',  so sigma_R includes the shared sigma_r only once.
	r, r_err = read_hist(inFile, f'hRatio_r_{q+1}_{x+1};1')
	if r is not None:
		pip = pip - r
		pim = pim - r

	# Mask invalid bins
	mask = (pip > 0) & (pim > 0)
	pip     = np.where(mask, pip,     np.nan)
	pim     = np.where(mask, pim,     np.nan)
	pip_err = np.where(mask, pip_err, np.nan)
	pim_err = np.where(mask, pim_err, np.nan)

	R = pip / pim
	if r is not None:
		r_err = np.where(mask, r_err, np.nan)
		R_err = (1 / pim) * np.sqrt(pip_err**2 + (R * pim_err)**2 + ((R - 1) * r_err)**2)
	else:
		R_err = R * np.sqrt((pip_err / pip)**2 + (pim_err / pim)**2)

	FF, FF_err = ratio_to_ff(R, R_err)

	return FF, FF_err, binCenters

#python makeRatios.py outfile file_1 title_1 file_2 title_2 ...
outFileName = sys.argv[1]
nFiles = int((len(sys.argv) - 2)/2)

inFile_list = [sys.argv[2*i + 2] for i in range(nFiles)]
inFile_tits = [sys.argv[2*i + 3] for i in range(nFiles)]


colorList = ['red', 'blue', 'green', 'magenta', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

for q in range(12):

	for x in range(14):

		fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')
		ax[0].set_ylim( [0, 1] )
		ax[1].set_xlim( [.3, .8] )
		ax[0].tick_params(axis='both', which='major', labelsize=14)
		ax[1].tick_params(axis='both', which='major', labelsize=14)

		z = np.linspace( .3, 1, 500 )
		ax[0].plot( z, np.asarray(ff(z)), color = 'black', linestyle = '--')
		ax[1].axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')

		ax[1].set_xlabel(r'$z$', fontsize=18)
		ax[1].set_ylabel(rf'$r(z)/r$({inFile_tits[0]})', fontsize=16)
		ax[0].set_ylabel(r'$r(z)$', fontsize=18)

		ax[0].text(.7, .75, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*(q), 2+.5*(q+1)), fontsize=18)
		ax[0].text(.7, .7, r'%.2f $< x_B <$ %.2f'%( .1 + .04*x, .1 + .04*(x+1)), fontsize=18)

		inFile = uproot.open(inFile_list[0])

		result = read_xs(inFile, q, x)
		if result is None:
			plt.close()
			continue

		values, errors, binCenters = result
		binCenters = binCenters + (-.01 + .002*q)

		if np.isnan(values).all():
			plt.close()
			continue

		ax[0].errorbar( binCenters, values, errors, label=inFile_tits[0], marker='.',
						color = colorList[0], linestyle = '',capsize = 2, lw = 1, capthick = 1, markersize=10,mec='black')

		rat_max = 1.05
		rat_min = 0.95
		for num in range( 1, len(inFile_list) ):
			inFile_temp = uproot.open(inFile_list[num])
			result_temp = read_xs(inFile_temp, q, x)
			if result_temp is None:
				continue

			values_temp, errors_temp, _ = result_temp

			ratio = values_temp / values
			ratio_err = ratio * np.sqrt( (errors / values)**2 + (errors_temp / values_temp)**2 )
			ratio_err[ratio_err <= 0] = np.nan

			ax[0].errorbar( binCenters, values_temp, errors_temp, label=inFile_tits[num], marker='.',
							color = colorList[num], linestyle = '',capsize = 2, lw = 1, capthick = 1,markersize=10,mec='black')

			ax[1].errorbar( binCenters, ratio, ratio_err, marker='.',
							color = colorList[num], linestyle = '',capsize = 2, lw = 1, capthick = 1, markersize=10,mec='black')

			if np.nanmax(ratio) > rat_max:
				rat_max = np.nanmax(ratio)*1.2
			if np.nanmin(ratio) < rat_min:
				rat_min = np.nanmin(ratio)*0.8

		ax[1].set_ylim( rat_min, rat_max )
		ax[0].legend()
		fig.savefig('{0}/ratio_{1}_{2}.png'.format(outFileName, q+1, x+1))

		plt.close()
