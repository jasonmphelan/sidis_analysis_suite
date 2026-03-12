import uproot
import numpy as np
import matplotlib.pyplot as plt
import sys

# ── Configuration ──────────────────────────────────────────────────────────────
# applyCorr value used when running makeRatioBinned / makeRatioBinned3D
#   > 3 : kaon correction applied (hRatio_k histograms present)
#   > 4 : rho  correction applied (hRatio_r histograms present)
APPLY_CORR = 4

# Number of extra variable bins (1 = no extra variable dimension).
# Overridden automatically when the file contains more bins.
BINS_VAR = 1

# True  = sum all var bins into a single curve (recommended for comparison plots)
# False = use var bin 0 only
SUM_VAR = True

# Kinematic bin definitions (must match cut_values.h)
BINS_Q2 = 12
BINS_XB = 14
Q2_MIN, Q2_MAX = 2.0, 8.0
XB_MIN, XB_MAX = 0.1, 0.66

# ── Helpers ────────────────────────────────────────────────────────────────────
def ff(z):
	return (1 - z) / (1 + z)

def get_hist(f, name):
	"""Return (values, errors, bin_centers) or None if histogram absent."""
	key = name + ";1"
	if key not in f:
		if name not in f:
			return None
		key = name
	h = f[key]
	edges = h.axis().edges()
	centers = (edges[:-1] + edges[1:]) / 2
	return np.array(h.values()), np.array(h.errors()), centers

def ratio_to_ff(R, R_err):
	"""Convert charge ratio R = pip/pim to asymmetry A = (4-R)/(4R-1)."""
	denom = 4 * R - 1
	with np.errstate(invalid='ignore', divide='ignore'):
		FF     = np.where(np.abs(denom) > 1e-9, (4 - R) / denom, np.nan)
		FF_err = np.where(np.abs(denom) > 1e-9, np.abs(15 / denom**2) * R_err, np.nan)
	return FF, FF_err

# ── Naming convention detection ────────────────────────────────────────────────
def detect_naming(inFile):
	"""Auto-detect histogram naming style and number of var bins in the file."""
	def _probe(name):
		return (name + ";1") in inFile or name in inFile
	def _probe_any(pattern_fn):
		for q in range(1, BINS_Q2 + 1):
			for x in range(1, BINS_XB + 1):
				if _probe(pattern_fn(q, x)):
					return True
		return False

	if _probe_any(lambda q, x: f"hRatio_0_{q}_{x}_1"):
		style = "new_multi"
		n_var = 0
		while _probe_any(lambda q, x: f"hRatio_{n_var}_{q}_{x}_1"):
			n_var += 1
	elif _probe_any(lambda q, x: f"hRatio_0_{q}_{x}"):
		style = "new_single"
		n_var = BINS_VAR
	elif _probe_any(lambda q, x: f"hRatio_{q}_{x}"):
		style = "old"
		n_var = 1
	else:
		keys = [k.split(";")[0] for k in inFile.keys() if "hRatio" in k][:8]
		raise RuntimeError(f"Cannot detect histogram naming convention. Sample keys: {keys}")
	return style, n_var

def make_hist_name(style):
	"""Return a hist_name(prefix, var, charge, q_idx, x_idx) -> str function."""
	def hist_name(prefix, var, charge, q_idx, x_idx):
		if style == "new_multi":
			charge_k = 0 if charge == '' else 1
			return prefix + f"_{var}{charge}_{q_idx+1}_{x_idx+1}_{charge_k+1}"
		elif style == "new_single":
			return prefix + f"_{var}{charge}_{q_idx+1}_{x_idx+1}"
		else:  # old
			pim_str = 'Pim' if charge == '_Pim' else ''
			return f"hRatio{pim_str}_{q_idx+1}_{x_idx+1}" if prefix == 'hRatio' \
				else prefix + f"{charge}_{q_idx+1}_{x_idx+1}"
	return hist_name

# ── Load all cross sections from one ROOT file ─────────────────────────────────
def load_file(path):
	"""
	Open a ROOT file, auto-detect histogram naming, load and correct pip/pim
	yields, apply kaon/rho corrections, sum over var bins (if configured),
	then compute A(z) = (4-R)/(4R-1).

	Returns
	-------
	FF     : ndarray, shape (BINS_Q2, BINS_XB, n_z)
	FF_err : ndarray, same shape
	z_cen  : ndarray, shape (n_z,)
	"""
	inFile = uproot.open(path)
	style, n_var = detect_naming(inFile)
	hist_name = make_hist_name(style)

	z_cen     = None
	n_z       = None
	pip_v_all = None
	pip_e_all = None
	pim_v_all = None
	pim_e_all = None

	for var in range(n_var):
		for q in range(BINS_Q2):
			for x in range(BINS_XB):
				pip = get_hist(inFile, hist_name('hRatio', var, '',     q, x))
				pim = get_hist(inFile, hist_name('hRatio', var, '_Pim', q, x))
				if pip is None or pim is None:
					continue

				pv, pe, cen = pip
				pmv, pme, _ = pim

				if z_cen is None:
					z_cen = cen
					n_z   = len(cen)
					shape = (n_var, BINS_Q2, BINS_XB, n_z)
					pip_v_all = np.full(shape, np.nan)
					pip_e_all = np.full(shape, np.nan)
					pim_v_all = np.full(shape, np.nan)
					pim_e_all = np.full(shape, np.nan)

				# Kaon correction (added to pip and pim independently)
				if APPLY_CORR > 3:
					k = get_hist(inFile, hist_name('hRatio_k', var, '',     q, x))
					if k is not None:
						pv  += k[0];  pe  = np.sqrt(pe**2  + k[1]**2)
					k = get_hist(inFile, hist_name('hRatio_k', var, '_Pim', q, x))
					if k is not None:
						pmv += k[0];  pme = np.sqrt(pme**2 + k[1]**2)

				# Rho correction (subtracted from pip and pim independently)
				if APPLY_CORR > 4:
					r = get_hist(inFile, hist_name('hRatio_r', var, '',     q, x))
					if r is not None:
						pv  -= r[0];  pe  = np.sqrt(pe**2  + r[1]**2)
					r = get_hist(inFile, hist_name('hRatio_r', var, '_Pim', q, x))
					if r is not None:
						pmv -= r[0];  pme = np.sqrt(pme**2 + r[1]**2)

				pip_v_all[var, q, x] = pv
				pip_e_all[var, q, x] = pe
				pim_v_all[var, q, x] = pmv
				pim_e_all[var, q, x] = pme

	if z_cen is None:
		return None, None, None

	# Sum over var bins (at yield level, before ratio)
	if n_var > 1 and SUM_VAR:
		any_pip = np.any(np.isfinite(pip_v_all), axis=0, keepdims=True)
		any_pim = np.any(np.isfinite(pim_v_all), axis=0, keepdims=True)
		pip_v_all = np.where(any_pip, np.nansum(pip_v_all, axis=0, keepdims=True), np.nan)
		pip_e_all = np.where(any_pip, np.sqrt(np.nansum(pip_e_all**2, axis=0, keepdims=True)), np.nan)
		pim_v_all = np.where(any_pim, np.nansum(pim_v_all, axis=0, keepdims=True), np.nan)
		pim_e_all = np.where(any_pim, np.sqrt(np.nansum(pim_e_all**2, axis=0, keepdims=True)), np.nan)

	# Collapse var axis to (BINS_Q2, BINS_XB, n_z)
	pip_v = pip_v_all[0]
	pip_e = pip_e_all[0]
	pim_v = pim_v_all[0]
	pim_e = pim_e_all[0]

	mask  = (pip_v > 0) & (pim_v > 0)
	pip_v = np.where(mask, pip_v, np.nan)
	pim_v = np.where(mask, pim_v, np.nan)
	pip_e = np.where(mask, pip_e, np.nan)
	pim_e = np.where(mask, pim_e, np.nan)

	R     = pip_v / pim_v
	R_err = R * np.sqrt((pip_e / pip_v)**2 + (pim_e / pim_v)**2)
	FF, FF_err = ratio_to_ff(R, R_err)

	return FF, FF_err, z_cen  # each shape (BINS_Q2, BINS_XB, n_z)

# ── Command-line interface ─────────────────────────────────────────────────────
# Usage: python plotRatios.py outdir file1 title1 [file2 title2 ...]
outFileName = sys.argv[1]
nFiles      = (len(sys.argv) - 2) // 2
inFile_list = [sys.argv[2*i + 2] for i in range(nFiles)]
inFile_tits = [sys.argv[2*i + 3] for i in range(nFiles)]

colorList = ['red', 'blue', 'green', 'magenta', 'brown', 'gold',
             'cyan', 'blueviolet', 'darkorange', 'black', 'yellow', 'gray'] * 3

# Pre-load all files up front
loaded = [load_file(p) for p in inFile_list]

z_line = np.linspace(0.3, 1.0, 500)

for q in range(BINS_Q2):
	for x in range(BINS_XB):

		FF0, FF0_err, z_cen = loaded[0]
		if FF0 is None or np.isnan(FF0[q, x]).all():
			print('none')
			continue

		fig, ax = plt.subplots(2, 1, figsize=(12, 6), height_ratios=[3, 1],
		                       sharex='col', layout='constrained')
		ax[0].set_ylim([0, 1])
		ax[1].set_xlim([0.3, 0.8])
		ax[0].tick_params(axis='both', which='major', labelsize=14)
		ax[1].tick_params(axis='both', which='major', labelsize=14)

		ax[0].plot(z_line, ff(z_line), color='black', linestyle='--')
		ax[1].axhline(y=1, color='black', linestyle='--')

		ax[1].set_xlabel(r'$z$', fontsize=18)
		ax[1].set_ylabel(rf'$r(z)/r$({inFile_tits[0]})', fontsize=16)
		ax[0].set_ylabel(r'$r(z)$', fontsize=18)

		q2_lo = Q2_MIN + q * (Q2_MAX - Q2_MIN) / BINS_Q2
		q2_hi = Q2_MIN + (q + 1) * (Q2_MAX - Q2_MIN) / BINS_Q2
		xb_lo = XB_MIN + x * (XB_MAX - XB_MIN) / BINS_XB
		xb_hi = XB_MIN + (x + 1) * (XB_MAX - XB_MIN) / BINS_XB
		ax[0].text(0.7, 0.75, rf'${q2_lo:.1f} < Q^2 < {q2_hi:.1f}$', fontsize=18)
		ax[0].text(0.7, 0.70, rf'${xb_lo:.2f} < x_B < {xb_hi:.2f}$', fontsize=18)

		vals0  = FF0[q, x]
		errs0  = FF0_err[q, x]
		z_plot = z_cen + (-0.01 + 0.002 * q)

		ax[0].errorbar(z_plot, vals0, errs0, label=inFile_tits[0], marker='.',
		               color=colorList[0], linestyle='', capsize=2, lw=1,
		               capthick=1, markersize=10, mec='black')

		rat_max, rat_min = 1.05, 0.95
		for num in range(1, nFiles):
			FF_n, FF_n_err, _ = loaded[num]
			if FF_n is None:
				continue
			vals_n = FF_n[q, x]
			errs_n = FF_n_err[q, x]

			ax[0].errorbar(z_plot, vals_n, errs_n, label=inFile_tits[num], marker='.',
			               color=colorList[num], linestyle='', capsize=2, lw=1,
			               capthick=1, markersize=10, mec='black')

			ratio     = vals_n / vals0
			ratio_err = ratio * np.sqrt((errs0 / vals0)**2 + (errs_n / vals_n)**2)
			ratio_err[ratio_err <= 0] = np.nan

			ax[1].errorbar(z_plot, ratio, ratio_err, marker='.',
			               color=colorList[num], linestyle='', capsize=2, lw=1,
			               capthick=1, markersize=10, mec='black')

			valid = ratio[np.isfinite(ratio)]
			if len(valid):
				rat_max = max(rat_max, np.nanmax(valid) * 1.2)
				rat_min = min(rat_min, np.nanmin(valid) * 0.8)

		ax[1].set_ylim(rat_min, rat_max)
		ax[0].legend()
		fig.savefig(f'{outFileName}/ratio_{q+1}_{x+1}.png')
		plt.close()
