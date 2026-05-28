#ifndef __CORRECTIONTOOLS__
#define __CORRECTIONTOOLS__

#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "TString.h"
#include "TH3F.h"
#include "TF1.h"
#include "TF3.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "constants.h"
#include "cut_values.h"
#define CORR_PATH _DATA

using namespace constants;
using namespace cutVals;

class correctionTools{
public:
	correctionTools(int _mode);
	~correctionTools();

	void setCorrectionLevel( int _mode ){ mode = _mode; }
	void loadParameters();
	void loadHistograms();

	void setWeightName( TString newName ){ weight_name = newName; };
	void setWeightName4D( TString newName ){ weight_name_4D = newName; };
	void setPi2kName( TString newName ){ pi_To_kaon_name = newName; };
	void setK2piName( TString newName ){ kaon_To_pi_name = newName; };

	// Continuous (TF3/TF1) fit files — set before calling loadContinuousCorrections()
	void setWeightFitName( TString newName ){ weight_fit_name = newName; }
	TString getWeightFitName() const { return weight_fit_name; }
	void setMisidFitName(  TString newName ){ misid_fit_name  = newName; }

	void setKinematics( double x, double q, double z, double p );
	void set4dBin(int inBin){ bin_4d = inBin; }
	void setN4dBins(int nBins ) { nBins_4d = nBins; }
	double getCorrectionFactor( int type, int charge );
	double getCorrectionError( int type, int charge );
//	void setFilePaths(int corr, TString path);
	void printFilePaths();
	void testHists(){ pi_to_k_Correction[0]->Print("all");}
	void loadFits();

	// Load continuous TF3 (mc correction) and TF1 (misid) fits.
	// When loaded, getCorrectionFactor/getCorrectionError use them in place of
	// the bin-centred TH3F/TH1F lookups for types 2, 3, and 4.
	void loadContinuousCorrections();

	void setCorrectionFile(TString name){ weight_name = name; }
	void loadNewEnergy( double energy );
	// Probe 4D correction file: returns bin count and populates corr_var_min/max/name
	int probeN4dBins();
	double getCorrVarMin()  { return corr_var_min; }
	double getCorrVarMax()  { return corr_var_max; }
	TString getCorrVarName(){ return corr_var_name; }

	double getX(){ return kin[0];}
	double getQ2(){ return kin[1];}
	double getZ(){ return kin[2];}
	double getP(){ return kin[3];}

	// ─────────────────────────────────────────────────────────────────────────
	// Full-uncertainty extensions (used by makeRatioBinned3D_fullUncertainty).
	//
	// These accessors expose the log-space MC-correction fit's
	//   - parameter vector  theta_mc^h    (size M)
	//   - covariance matrix V_theta_mc^h  (M x M)
	//   - basis vector      f_mc^h(xB, Q2, z)  (size M)
	//
	// where ln C_mc^h(xB, Q2, z) = theta_mc^h . f_mc^h(xB, Q2, z)
	// and the basis ordering is (block-major):
	//
	//   Block A — polynomial in (ln xB, ln Q2):
	//     for s = 0..dP, for i = 0..s, j = s - i :  (ln xB)^i * (ln Q2)^j
	//     size = (dP+1)*(dP+2)/2
	//
	//   Block B — z * polynomial in ln Q2:
	//     for k = 1..dR, for m = 0..dZQ :  (ln Q2)^m * z^k
	//     size = dR * (dZQ + 1)
	//
	//   Block C — z * polynomial in ln xB:
	//     for k = 1..dR, for n = 1..dZX :  (ln xB)^n * z^k
	//     size = dR * dZX
	//
	// Total M = sizeA + sizeB + sizeC.
	//
	// Storage convention in the corrections fit ROOT file (loaded by
	// loadContinuousCorrections / loadMcFitCovariances):
	//   theta_McCorrectionP / theta_McCorrectionM   (TVectorD of length M)
	//   cov_McCorrectionP   / cov_McCorrectionM     (TMatrixDSym, M x M)
	//   dP_McCorrection, dR_McCorrection,
	//   dZQ_McCorrection, dZX_McCorrection          (TParameter<int>)
	//
	// If any of these objects is missing, the corresponding accessor
	// returns false and the caller must fall back to the legacy
	// scalar-sigma treatment.
	// ─────────────────────────────────────────────────────────────────────────
	void loadMcFitCovariances();          // called automatically by loadContinuousCorrections
	bool hasMcFitCovariances(int charge) const;
	int  getMcNParams(int charge) const;
	bool getMcParams(int charge, std::vector<double>& theta_out) const;
	bool getMcCov   (int charge, TMatrixDSym&        V_out)     const;
	bool evaluateMcBasis(double xB, double Q2, double z,
	                     int charge, std::vector<double>& f_out) const;
	int  getMcBasisDeg_dP () const { return mc_dP;  }
	int  getMcBasisDeg_dR () const { return mc_dR;  }
	int  getMcBasisDeg_dZQ() const { return mc_dZQ; }
	int  getMcBasisDeg_dZX() const { return mc_dZX; }

	// Static helper: fills f with the basis vector at (xB, Q2, z) for the
	// given polynomial degrees, using the ordering documented above.
	static void buildMcBasis(double xB, double Q2, double z,
	                         int dP, int dR, int dZQ, int dZX,
	                         std::vector<double>& f);
	static int  mcBasisSize(int dP, int dR, int dZQ, int dZX){
		return (dP+1)*(dP+2)/2 + dR*(dZQ+1) + dR*dZX;
	}



private:
	int 	mode;
	TString weight_name     = "corrections_10.2_AN_test.root";
	TString weight_name_4D  = "corrections_10.2_phi_q.root";
	TString kaon_To_pi_name = "corrections_misid_AN.root";
	TString pi_To_kaon_name = "corrections_misid_AN.root";

	// Continuous fit files (written by plotCorrections3D.ipynb / plotMisidCorrections.ipynb)
	TString weight_fit_name = "corrections_fit.root";
	TString misid_fit_name  = "corrections_misid_fit.root";

	double kin[4];
	int bin_4d = -1;
	int nBins_4d = 0;
	double misid_min = 2;
	double corr_var_min  = 0;
	double corr_var_max  = 1;
	TString corr_var_name = "";

	// ── Binned histograms (original) ──────────────────────────────────────────
	TH3F *	accCorrection[2]      = {nullptr, nullptr};
	TH3F *	binMigration[2]       = {nullptr, nullptr};
	TH3F *	mcCorrection[2]       = {nullptr, nullptr};
	TH1F *	k_to_pi_Correction[2] = {nullptr, nullptr};
	TH1F *	pi_to_k_Correction[2] = {nullptr, nullptr};

	// ── Continuous fits (loaded by loadContinuousCorrections) ─────────────────
	// MC correction:  TF3(xB, Q², z) — hMcCorrectionP/M_tf3
	// Sigma on MC:    TF3(xB, Q², z) — hMcCorrectionP/M_tf3_sigma
	// 4D MC corr:     TF3(xB, Q², z) — hMcCorrectionP/M_tf3_{bin} (one per 4th-var bin)
	// Misid ratio:    TF1(p)         — hRatio_p_{pip/pim}_tf1  (pi→k)
	//                                  hRatio_p_{kp/km}_tf1    (k→pi)
	// Sigma on misid: TF1(p)         — corresponding _tf1_sigma
	TF3 *	mcCorrectionFit[2]        = {nullptr, nullptr};
	TF3 *	mcCorrectionFitSigma[2]   = {nullptr, nullptr};
	TF3 *	mcCorrectionFit4D[2][10]  = {};
	TF3 *	mcCorrectionFit4DSigma[2][10] = {};
	TF1 *	pi2kRatioFit[2]         = {nullptr, nullptr};
	TF1 *	pi2kRatioFitSigma[2]    = {nullptr, nullptr};
	TF1 *	k2piRatioFit[2]         = {nullptr, nullptr};
	TF1 *	k2piRatioFitSigma[2]    = {nullptr, nullptr};

	// ── Log-space-fit theta and covariance (full-uncertainty path) ────────────
	// mc_theta[charge]  : parameter vector of length M
	// mc_cov[charge]    : M x M parameter covariance
	// mc_dP/dR/dZQ/dZX  : basis polynomial degrees (see comment block above)
	// mc_cov_loaded[c]  : true iff theta + cov + degrees were all loaded
	std::vector<double> mc_theta[2];
	TMatrixDSym         mc_cov  [2];
	int  mc_dP  = -1;
	int  mc_dR  = -1;
	int  mc_dZQ = -1;
	int  mc_dZX = -1;
	bool mc_cov_loaded[2] = {false, false};

	// ── Legacy per-bin fits ───────────────────────────────────────────────────
	TF1 * weightFit[2][bins_xB][bins_Q2];
	TF1 * pi2kFit[2][bins_xB][bins_Q2][4];
	TF1 * k2piFit[2][bins_xB][bins_Q2][4];

	double weight_Parameters[2][4][4][4];
	double k_To_pi_Parameters[2][4][5][4][4];
	double pi_To_k_Parameters[2][4][5][4][4];

	void loadWeightParameters();
	void loadWeightParameters_4D();
	void loadkTopiParameters();
	void loadpiTokParameters();

};

#endif
