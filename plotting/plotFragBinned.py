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

def _probe(f, name):
	return (name + ";1") in f or name in f

def detect_style(f):
	"""Return 'ratio' if hRatio_0_q_x exists (new code with covariance-aware errors),
	'sigma_single' if only hSigmaPlus/hSigmaMinus_Pim are present, else 'ratio'."""
	for q in range(1, 13):
		for x in range(1, 15):
			if _probe(f, f'hRatio_0_{q}_{x}'):
				return 'ratio'
			if _probe(f, f'hSigmaPlus_0_{q}_{x}'):
				return 'sigma_single'
	return 'ratio'

def load_ratio(f, style, q, x):
	"""Return (values, errors, xEdges) for bin (q+1, x+1), or None if empty."""
	if style == 'sigma_single':
		kp = f'hSigmaPlus_0_{q+1}_{x+1}'
		km = f'hSigmaMinus_Pim_0_{q+1}_{x+1}'
		if not (_probe(f, kp) and _probe(f, km)):
			return None
		hp = f[kp]; hm = f[km]
		pv = np.array(hp.values()); pe = np.array(hp.errors())
		mv = np.array(hm.values()); me = np.array(hm.errors())
		mask = (pv > 0) & (mv > 0)
		pv = np.where(mask, pv, np.nan); mv = np.where(mask, mv, np.nan)
		R = pv / mv
		R_err = R * np.sqrt((pe / pv)**2 + (me / mv)**2)
		return R, R_err, hp.axis().edges()
	else:
		key = f'hRatio_0_{q+1}_{x+1}'
		if not _probe(f, key):
			return None
		h = f[key]
		return np.array(h.values()), np.array(h.errors()), h.axis().edges()

#python plotFragBinned.py outfile infile [mode]
# mode: 'x' = loop over xB bins, color by Q2 (default)
#        'q' = loop over Q2 bins, color by xB
outFileName = sys.argv[1]
inFile = uproot.open(sys.argv[2])
mode = sys.argv[3] if len(sys.argv) > 3 else 'x'
style = detect_style(inFile)

colorList = ['red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray', 'red', 'blue', 'magenta', 'green', 'brown', 'gold', 'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray']

if mode == 'q':
	outer_range = range(12)
else:
	outer_range = range(14)

for outer in outer_range:

	fig, ax = plt.subplots(2, 1, figsize=(12,6), height_ratios=[3,1], sharex='col', layout='constrained')

	ax[0].set_ylim( [0, 5] )
	ax[1].set_xlim( [.3, .8] )

	ax[0].tick_params(axis='both', which='major', labelsize=14)
	ax[1].tick_params(axis='both', which='major', labelsize=14)

	z = np.linspace( .3, 1, 500 )
	ax[0].plot( z, np.asarray(ff(z)), color = 'black', linestyle = '--')
	ax[1].axhline( y=1, xmin=0, xmax=1 , color = 'black', linestyle = '--')

	ax[1].set_xlabel(r'$z$', fontsize=18)
	ax[1].set_ylabel(rf'$r(z, Q^2)/$FF(z))', fontsize=18)
	ax[0].set_ylabel(r'$r(z)$', fontsize=18)

	if mode == 'q':
		ax[0].text(.7, .7, r'%.1f $< Q^2 <$ %.1f'%(2 + .5*outer, 2 + .5*(outer+1)), fontsize=18)
		inner_range = range(14)
	else:
		ax[0].text(.7, .7, r'%.2f $< x_B <$ %.2f'%(.1 + .04*outer, .1 + .04*(outer+1)), fontsize=18)
		inner_range = range(12)

	rat_max = 1.05
	rat_min = 0.95
	first_bin = -1
	for inner in inner_range:

		if mode == 'q':
			q, x = outer, inner
		else:
			q, x = inner, outer

		result = load_ratio(inFile, style, q, x)
		if result is None:
			continue
		R, R_err, xEdges = result

		R[R <= 0] = np.nan
		if np.isnan(R).all():
			continue
		if first_bin < 0:
			first_bin = inner

		denom = 4 * R - 1
		values = R#np.where(np.abs(denom) > 1e-9, (4 - R) / denom, np.nan)
		errors = R_err#np.where(np.abs(denom) > 1e-9, np.abs(15 / denom**2) * R_err, np.nan)
		binCenters = np.asarray((xEdges[:-1]+xEdges[1:])/2 + (.002*(inner-first_bin)))

		if mode == 'q':
			label = r'%.2f $< x_B <$ %.2f'%(.1 + .04*x, .1 + .04*(x+1))
		else:
			label = r'%.1f $< Q^2 <$ %.1f'%(2 + .5*q, 2 + .5*(q+1))

		ax[0].errorbar(binCenters, values, errors, label=label, marker='o',
					color=colorList[inner], linestyle='', capsize=2, lw=1, capthick=1, mec='black')

		ratio = values/np.asarray(ff(binCenters))
		ratio_err = ratio*np.sqrt((errors/values)**2)
		ratio_err[ratio_err<=0] = np.nan

		ax[1].errorbar(binCenters, ratio, ratio_err, marker='o',
					color=colorList[inner], linestyle='', capsize=2, lw=1, capthick=1, mec='black')

		if np.nanmax(ratio) > rat_max:
			rat_max = np.nanmax(ratio)*1.2
		if np.nanmin(ratio) < rat_min:
			rat_min = np.nanmin(ratio)*0.8

		for edge in xEdges:
			ax[0].axvline(x=edge, color='gray', linestyle=':', lw=0.8, alpha=0.5)
			ax[1].axvline(x=edge, color='gray', linestyle=':', lw=0.8, alpha=0.5)

	ax[1].set_ylim([0, 3])
	ax[0].legend()
	fig.savefig('{0}/ratio_{1}.png'.format(outFileName, outer+1))

	plt.close()
