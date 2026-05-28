import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description="Plot pi+/pi- charge comparison across ROOT files.")
parser.add_argument("files", nargs='+', help="Input ROOT files")
parser.add_argument("-o", "--output", required=True, help="Output directory")
parser.add_argument("-l", "--labels", nargs='+', help="Legend labels (one per file, defaults to filename stems)")
args = parser.parse_args()

if args.labels and len(args.labels) != len(args.files):
	parser.error(f"Number of labels ({len(args.labels)}) must match number of files ({len(args.files)})")

out_dir = args.output
labels = args.labels if args.labels else [Path(f).stem for f in args.files]
inFiles = [(label, uproot.open(f)) for label, f in zip(labels, args.files)]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def setTitle2D(ax, key):
	if "Beta" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]")
		ax.set_ylabel(r"$\beta$")
	if "WV" in key:
		ax.set_xlabel(r"V$_{PCAL}$ [cm]")
		ax.set_ylabel(r"W$_{PCAL}$ [cm]")
	if "Edep" in key:
		ax.set_xlabel(r"$E_{PCAL}$ [GeV]")
		ax.set_ylabel(r"E_{CIN} + E_{COUT} [GeV]")
	if "SF_corr" in key:
		ax.set_xlabel(r"$\delta E_{dep}(PCAL)/p_{e}$")
		ax.set_ylabel(r"$(E_{PCAL} + E_{CIN} + E_{COUT})/p_e$ [GeV]")
	if "SF_sec" in key:
		ax.set_xlabel(r"$p_{e}$ [GeV]")
		ax.set_ylabel(r"$\delta E_{CIN}/p_e$")
	if "Chi2" in key:
		ax.set_xlabel(r"$p_{\pi}$ [GeV]")
		ax.set_ylabel(r"$\chi^2_{PID}$")

def setXLabel(ax, key):
	labels = {
		"Beta":    r"$\beta$",
		"Vz_e":    r"$v_z^e$ [cm]",
		"Vz_pi":   r"$v_z^{\pi} - v_z^e$ [cm]",
		"Z_pi":    r"$z$",
		"Y_pi":    r"$y$",
		"W_pi":    r"$W$ [GeV]",
		"Theta_e": r"$\theta_e$ [deg.]",
		"Theta_pi":r"$\theta_{\pi}$ [deg.]",
		"Q2":      r"$Q^2$ [GeV$^2$]",
		"P_pi":    r"$p_{\pi}$ [GeV]",
		"P_e":     r"$p_{e}$ [GeV]",
		"Mx":      r"$M_X$ [GeV]",
		"Omega":   r"$\omega$ [GeV]",
		"Xb":      r"$x_B$",
		"Xf":      r"$x_F$",
		"Eta":     r"$\eta$",
		"Pt":      r"$p_{\pi}^{\perp}$ [GeV]",
		"Phi_q":   r"$\phi_{\pi}^{q\text{-frame}}$ [deg.]",
	}
	for substr, label in labels.items():
		if substr in key:
			ax.set_xlabel(label, fontsize=18)
			break

# Drive loop from keys in the first file, skipping pim keys
ref_label, ref_file = inFiles[0]

for key in ref_file.keys():
	if "pim" in key:
		continue
	if "0_0" not in key:
		continue
	hist = ref_file[key]

	if not isinstance(hist, uproot.models.TH.Model_TH1F_v3):
		continue
	
	hName_pip = key
	hName_pim = key.replace("pip", "pim")

	fig, axes = plt.subplots(2, 1, figsize=(10, 8),
		height_ratios=[3, 2], sharex='col', layout='constrained')

	for i, (label, f) in enumerate(inFiles):
		color_pip = colors[(2 * i) % len(colors)]
		color_pim = colors[(2 * i + 1) % len(colors)]

		pip = f[hName_pip]
		pim = f[hName_pim]

		values_pip = pip.values().copy()
		errors_pip = pip.errors().copy()
		values_pim = pim.values().copy()
		errors_pim = pim.errors().copy()

		norm_pip = np.nansum(values_pip)
		norm_pim = np.nansum(values_pim)

		values_pip[values_pip == 0] = np.nan
		values_pim[values_pim == 0] = np.nan

		ratio = values_pip / values_pim
		ratio_err = ratio * np.sqrt((errors_pip / values_pip)**2 + (errors_pim / values_pim)**2)

		values_pip = values_pip / norm_pip
		errors_pip = errors_pip / norm_pip
		values_pim = values_pim / norm_pim
		errors_pim = errors_pim / norm_pim

		xEdges = pip.axis().edges()
		binCenters = (xEdges[:-1] + xEdges[1:]) / 2

		axes[0].errorbar(binCenters, values_pip, errors_pip, capsize=2,
			fmt='-', color=color_pip, ecolor=color_pip, drawstyle='steps-mid', label=r"$d(e, e'\pi^+$)" )
		axes[0].errorbar(binCenters, values_pim, errors_pim, capsize=2,
			fmt='-', color=color_pim, ecolor=color_pim, drawstyle='steps-mid', label=r"$d(e, e'\pi^-$)")
		axes[1].errorbar(binCenters, ratio, ratio_err, capsize=2,
			fmt='o', color=color_pip, ecolor=color_pip, label=label)

	axes[1].axhline(y=1, linestyle='--', color='grey')

	for ax in axes:
		ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
		ax.tick_params(which='both', width=2, labelsize=12)
		ax.tick_params(which='major', length=7)
		ax.tick_params(which='minor', length=4)

	axes[0].set_ylabel(r"Yield [a.u.]", fontsize=14)
	axes[1].set_ylabel(r"$\pi^+/\pi^-$", fontsize=14)
	setXLabel(axes[1], key)

	pdfName = hName_pip.replace("_pip", "")
	axes[0].legend(fontsize=12, loc='best')
	fig.savefig(f"{out_dir}/{pdfName}.png")
	plt.close(fig)
