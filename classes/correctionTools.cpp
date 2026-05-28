#include "correctionTools.h"
#include "TFile.h"
#include "TF1.h"
#include "TF3.h"
#include "TH3F.h"
#include "TParameter.h"
#include "TNamed.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include <cmath>
#include "constants.h"
#include "cut_values.h"
#define CORR_PATH _DATA

using namespace constants;
using namespace cutVals;

correctionTools::correctionTools(int _mode){
    setCorrectionLevel (_mode);
}

correctionTools::~correctionTools(){}
	
void correctionTools::printFilePaths(){
	std::cout<< (TString) _DATA + "/correctionFiles/"+ weight_name<<std::endl;
	std::cout<< (TString) _DATA + "/correctionFiles/"+ weight_name_4D<<std::endl;

	std::cout<< (TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name<<std::endl;
	std::cout<< (TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name<<std::endl;
}

void correctionTools::loadHistograms(){
//	TFile weightFile((TString) _DATA + "/correctionFiles/"+ weight_name);
//	TFile pi2kFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
//	TFile k2piFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	//weightFile->Close();
	//pi2kFile->Close();
	//k2piFile->Close();
	
	TFile * weightFile = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	TFile * weightFile_4d = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name_4D);

	TFile * pi2kFile = new TFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	TFile * k2piFile = new TFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	
	accCorrection[0] = (TH3F*)weightFile->Get( "hAccCorrectionP" );
	accCorrection[1] = (TH3F*)weightFile->Get( "hAccCorrectionM" );
	binMigration[0] = (TH3F*)weightFile->Get( "hBinMigrationP" );
	binMigration[1] = (TH3F*)weightFile->Get( "hBinMigrationM" );
	mcCorrection[0] = (TH3F*)weightFile->Get( "hMcCorrectionP" );
	mcCorrection[1] = (TH3F*)weightFile->Get( "hMcCorrectionM" );


	accCorrection[0]->SetDirectory(0);
	accCorrection[1]->SetDirectory(0);
	binMigration[0]->SetDirectory(0);
	binMigration[1]->SetDirectory(0);
	mcCorrection[0]->SetDirectory(0);
	mcCorrection[1]->SetDirectory(0);

	k_to_pi_Correction[0] = (TH1F*)pi2kFile->Get( "hFactor_kp" );
	k_to_pi_Correction[1] = (TH1F*)pi2kFile->Get( "hFactor_km" );
	pi_to_k_Correction[0] = (TH1F*)pi2kFile->Get( "hFactor_pip" );
	pi_to_k_Correction[1] = (TH1F*)pi2kFile->Get( "hFactor_pim" );

	if( !k_to_pi_Correction[0] || !k_to_pi_Correction[1] ||
	    !pi_to_k_Correction[0] || !pi_to_k_Correction[1] ){
		std::cerr << "ERROR: Missing kaon correction histograms in " << pi_To_kaon_name << "\n";
		std::cerr << "  hFactor_kp:  " << (k_to_pi_Correction[0] ? "OK" : "MISSING") << "\n";
		std::cerr << "  hFactor_km:  " << (k_to_pi_Correction[1] ? "OK" : "MISSING") << "\n";
		std::cerr << "  hFactor_pip: " << (pi_to_k_Correction[0] ? "OK" : "MISSING") << "\n";
		std::cerr << "  hFactor_pim: " << (pi_to_k_Correction[1] ? "OK" : "MISSING") << "\n";
	}

	if( k_to_pi_Correction[0] ) k_to_pi_Correction[0]->SetDirectory(0);
	if( k_to_pi_Correction[1] ) k_to_pi_Correction[1]->SetDirectory(0);
	if( pi_to_k_Correction[0] ) pi_to_k_Correction[0]->SetDirectory(0);
	if( pi_to_k_Correction[1] ) pi_to_k_Correction[1]->SetDirectory(0);

	
	/*
	for( int i = 0; i < 4; i++ ){
	
		k_to_pi_Correction[0][i] = (TH3F*)pi2kFile->Get( Form( "hKaonCorrP_%i", i) );
		k_to_pi_Correction[1][i] = (TH3F*)pi2kFile->Get( Form( "hKaonCorrM_%i", i) );
		pi_to_k_Correction[0][i] = (TH3F*)k2piFile->Get( Form( "hKaonCorrP_%i", i) );
		pi_to_k_Correction[1][i] = (TH3F*)k2piFile->Get( Form( "hKaonCorrM_%i", i) );

		k_to_pi_Correction[0][i]->SetDirectory(0);
		k_to_pi_Correction[1][i]->SetDirectory(0);
		pi_to_k_Correction[0][i]->SetDirectory(0);
		pi_to_k_Correction[1][i]->SetDirectory(0);
	}
	*/
	for( int i = 0; i < nBins_4d; i++ ){
		auto* tf3P = (TF3*)weightFile_4d->Get( Form( "hMcCorrectionP_tf3_%i", i ) );
		auto* tf3M = (TF3*)weightFile_4d->Get( Form( "hMcCorrectionM_tf3_%i", i ) );
		if( tf3P ) mcCorrectionFit4D[0][i] = (TF3*)tf3P->Clone();
		if( tf3M ) mcCorrectionFit4D[1][i] = (TF3*)tf3M->Clone();
	}

	//delete weightFile;
	//delete pi2kFile;
	//delete k2piFile;
	delete weightFile;
	delete weightFile_4d;

	delete pi2kFile;
	delete k2piFile;

}

int correctionTools::probeN4dBins(){
	TFile * f = new TFile((TString) _DATA + "/correctionFiles/" + weight_name_4D);
	if( !f || f->IsZombie() ){ delete f; return 0; }
	int count = 0;
	while( f->Get( Form("hMcCorrP_fit_%d", count) ) ) count++;
	if( count > 0 ){
		auto* pMin  = (TParameter<double>*)f->Get("var_min");
		auto* pMax  = (TParameter<double>*)f->Get("var_max");
		auto* pName = (TNamed*)f->Get("var_name");
		if( pMin  ) corr_var_min  = pMin->GetVal();
		if( pMax  ) corr_var_max  = pMax->GetVal();
		if( pName ) corr_var_name = pName->GetTitle();
	}
	delete f;
	return count;
}

void correctionTools::loadNewEnergy( double energy ){
	setWeightName(Form("corrections_%.1f_fit.root", energy));
	loadFits();
}

void correctionTools::loadContinuousCorrections(){
	TString weightFitPath = (TString)_DATA + "/correctionFiles/" + weight_fit_name;
	TString misidFitPath  = (TString)_DATA + "/correctionFiles/" + misid_fit_name;

	// ── MC correction TF3 ────────────────────────────────────────────────────
	TFile * wf = new TFile(weightFitPath);
	if( !wf || wf->IsZombie() ){
		std::cerr << "WARNING: could not open " << weightFitPath
		          << " — continuous MC corrections unavailable\n";
	} else {
		// Clone() gives ownership-independent heap copies; safe to delete wf after
		auto* tf3p  = (TF3*)wf->Get("hMcCorrectionP_tf3");
		auto* tf3m  = (TF3*)wf->Get("hMcCorrectionM_tf3");
		auto* tf3ps = (TF3*)wf->Get("hMcCorrectionP_tf3_sigma");
		auto* tf3ms = (TF3*)wf->Get("hMcCorrectionM_tf3_sigma");
		if( tf3p  ) mcCorrectionFit[0]      = (TF3*)tf3p ->Clone();
		if( tf3m  ) mcCorrectionFit[1]      = (TF3*)tf3m ->Clone();
		if( tf3ps ) mcCorrectionFitSigma[0] = (TF3*)tf3ps->Clone();
		if( tf3ms ) mcCorrectionFitSigma[1] = (TF3*)tf3ms->Clone();

		bool ok = mcCorrectionFit[0] && mcCorrectionFit[1];
		std::cout << "loadContinuousCorrections: MC TF3 "
		          << (ok ? "loaded OK" : "MISSING") << " from " << weightFitPath << "\n";

		// ── 4D TF3 corrections (one per 4th-variable bin) ───────────────────────
		{
			int count = 0;
			while( wf->Get( Form("hMcCorrP_fit_%d", count) ) ) count++;
			if( count > 0 ){
				nBins_4d = count;
				for( int i = 0; i < nBins_4d; i++ ){
					auto* tf3P  = (TF3*)wf->Get( Form("hMcCorrectionP_tf3_%d",       i) );
					auto* tf3M  = (TF3*)wf->Get( Form("hMcCorrectionM_tf3_%d",       i) );
					auto* tf3PS = (TF3*)wf->Get( Form("hMcCorrectionP_tf3_sigma_%d", i) );
					auto* tf3MS = (TF3*)wf->Get( Form("hMcCorrectionM_tf3_sigma_%d", i) );
					if( tf3P  ) mcCorrectionFit4D[0][i]      = (TF3*)tf3P ->Clone();
					if( tf3M  ) mcCorrectionFit4D[1][i]      = (TF3*)tf3M ->Clone();
					if( tf3PS ) mcCorrectionFit4DSigma[0][i] = (TF3*)tf3PS->Clone();
					if( tf3MS ) mcCorrectionFit4DSigma[1][i] = (TF3*)tf3MS->Clone();
				}
				bool ok = mcCorrectionFit4D[0][0] && mcCorrectionFit4D[1][0];
				std::cout << "loadContinuousCorrections: 4D TF3 bins "
				          << (ok ? "loaded OK" : "MISSING") << " ("
				          << nBins_4d << " bins) from " << weightFitPath << "\n";
			}
		}

		wf->Close();
		delete wf;
	}

	// ── Misid ratio TF1 ───────────────────────────────────────────────────────
	// pi→k  (type 3):  hRatio_p_pip_tf1 / hRatio_p_pim_tf1
	// k→pi  (type 4):  hRatio_p_kp_tf1  / hRatio_p_km_tf1
	TFile * mf = new TFile(misidFitPath);
	if( !mf || mf->IsZombie() ){
		std::cerr << "WARNING: could not open " << misidFitPath
		          << " — continuous misid corrections unavailable\n";
	} else {
		// Clone() gives ownership-independent heap copies; safe to delete mf after
		auto* f_pip  = (TF1*)mf->Get("hRatio_p_pip_tf1");
		auto* f_pim  = (TF1*)mf->Get("hRatio_p_pim_tf1");
		auto* f_pips = (TF1*)mf->Get("hRatio_p_pip_tf1_sigma");
		auto* f_pims = (TF1*)mf->Get("hRatio_p_pim_tf1_sigma");
		auto* f_kp   = (TF1*)mf->Get("hRatio_p_kp_tf1");
		auto* f_km   = (TF1*)mf->Get("hRatio_p_km_tf1");
		auto* f_kps  = (TF1*)mf->Get("hRatio_p_kp_tf1_sigma");
		auto* f_kms  = (TF1*)mf->Get("hRatio_p_km_tf1_sigma");
		if( f_pip  ) pi2kRatioFit[0]      = (TF1*)f_pip ->Clone();
		if( f_pim  ) pi2kRatioFit[1]      = (TF1*)f_pim ->Clone();
		if( f_pips ) pi2kRatioFitSigma[0] = (TF1*)f_pips->Clone();
		if( f_pims ) pi2kRatioFitSigma[1] = (TF1*)f_pims->Clone();
		if( f_kp   ) k2piRatioFit[0]      = (TF1*)f_kp  ->Clone();
		if( f_km   ) k2piRatioFit[1]      = (TF1*)f_km  ->Clone();
		if( f_kps  ) k2piRatioFitSigma[0] = (TF1*)f_kps ->Clone();
		if( f_kms  ) k2piRatioFitSigma[1] = (TF1*)f_kms ->Clone();

		bool ok = pi2kRatioFit[0] && pi2kRatioFit[1] &&
		          k2piRatioFit[0] && k2piRatioFit[1];
		std::cout << "loadContinuousCorrections: misid TF1 "
		          << (ok ? "loaded OK" : "MISSING") << " from " << misidFitPath << "\n";
		mf->Close();
		delete mf;
	}

	// Try to load the log-space fit's theta vector and covariance matrix.
	// This is optional: the legacy scalar-sigma path still works if these
	// objects are absent in the fit file.  It is required by the
	// makeRatioBinned3D_fullUncertainty executable.
	loadMcFitCovariances();
}

// ─────────────────────────────────────────────────────────────────────────────
// Full-uncertainty extensions
// ─────────────────────────────────────────────────────────────────────────────

void correctionTools::loadMcFitCovariances(){
	TString weightFitPath = (TString)_DATA + "/correctionFiles/" + weight_fit_name;
	TFile * wf = new TFile(weightFitPath);
	if( !wf || wf->IsZombie() ){
		std::cerr << "loadMcFitCovariances: could not open " << weightFitPath
		          << " — log-space-fit covariance unavailable\n";
		if( wf ) delete wf;
		return;
	}

	// Basis polynomial degrees — required to interpret theta / cov.
	auto* p_dP  = (TParameter<int>*) wf->Get("dP_McCorrection");
	auto* p_dR  = (TParameter<int>*) wf->Get("dR_McCorrection");
	auto* p_dZQ = (TParameter<int>*) wf->Get("dZQ_McCorrection");
	auto* p_dZX = (TParameter<int>*) wf->Get("dZX_McCorrection");
	if( p_dP && p_dR && p_dZQ && p_dZX ){
		mc_dP  = p_dP ->GetVal();
		mc_dR  = p_dR ->GetVal();
		mc_dZQ = p_dZQ->GetVal();
		mc_dZX = p_dZX->GetVal();
	} else {
		std::cerr << "loadMcFitCovariances: basis-degree TParameters missing in "
		          << weightFitPath
		          << " (expected dP_McCorrection, dR_McCorrection, "
		          << "dZQ_McCorrection, dZX_McCorrection)\n";
		wf->Close();
		delete wf;
		return;
	}
	const int M = mcBasisSize(mc_dP, mc_dR, mc_dZQ, mc_dZX);

	const char* theta_keys[2] = {"theta_McCorrectionP", "theta_McCorrectionM"};
	const char* cov_keys  [2] = {"cov_McCorrectionP",   "cov_McCorrectionM"};

	for( int c = 0; c < 2; c++ ){
		mc_cov_loaded[c] = false;
		auto* tv = (TVectorD*) wf->Get(theta_keys[c]);
		auto* tc = (TMatrixDSym*) wf->Get(cov_keys[c]);
		if( !tv || !tc ){
			std::cerr << "loadMcFitCovariances: missing "
			          << theta_keys[c] << " or " << cov_keys[c]
			          << " in " << weightFitPath << "\n";
			continue;
		}
		if( tv->GetNrows() != M || tc->GetNrows() != M || tc->GetNcols() != M ){
			std::cerr << "loadMcFitCovariances: dimension mismatch for charge "
			          << c << ": theta has " << tv->GetNrows()
			          << ", cov is " << tc->GetNrows() << "x" << tc->GetNcols()
			          << ", but basis size = " << M
			          << " (dP=" << mc_dP << ", dR=" << mc_dR
			          << ", dZQ=" << mc_dZQ << ", dZX=" << mc_dZX << ")\n";
			continue;
		}
		mc_theta[c].assign(M, 0.0);
		for( int a = 0; a < M; a++ ) mc_theta[c][a] = (*tv)[a];
		mc_cov[c].ResizeTo(M, M);
		mc_cov[c] = *tc;
		mc_cov_loaded[c] = true;
	}

	std::cout << "loadMcFitCovariances: theta+cov "
	          << (mc_cov_loaded[0] ? "OK" : "MISSING") << " (+), "
	          << (mc_cov_loaded[1] ? "OK" : "MISSING") << " (-), "
	          << "M=" << M
	          << " (dP=" << mc_dP << ", dR=" << mc_dR
	          << ", dZQ=" << mc_dZQ << ", dZX=" << mc_dZX << ")\n";

	wf->Close();
	delete wf;
}

bool correctionTools::hasMcFitCovariances(int charge) const{
	if( charge != 0 && charge != 1 ) return false;
	return mc_cov_loaded[charge];
}

int correctionTools::getMcNParams(int charge) const{
	if( charge != 0 && charge != 1 ) return 0;
	return mc_cov_loaded[charge] ? (int)mc_theta[charge].size() : 0;
}

bool correctionTools::getMcParams(int charge, std::vector<double>& theta_out) const{
	if( charge != 0 && charge != 1 ) return false;
	if( !mc_cov_loaded[charge] ) return false;
	theta_out = mc_theta[charge];
	return true;
}

bool correctionTools::getMcCov(int charge, TMatrixDSym& V_out) const{
	if( charge != 0 && charge != 1 ) return false;
	if( !mc_cov_loaded[charge] ) return false;
	V_out.ResizeTo(mc_cov[charge].GetNrows(), mc_cov[charge].GetNcols());
	V_out = mc_cov[charge];
	return true;
}

bool correctionTools::evaluateMcBasis(double xB, double Q2, double z,
                                      int charge, std::vector<double>& f_out) const{
	if( charge != 0 && charge != 1 ) return false;
	if( mc_dP < 0 || mc_dR < 0 || mc_dZQ < 0 || mc_dZX < 0 ) return false;
	buildMcBasis(xB, Q2, z, mc_dP, mc_dR, mc_dZQ, mc_dZX, f_out);
	return true;
}

// Static helper: writes the basis monomials into f, in the documented order.
//   Block A:  for s=0..dP, for i=0..s, j=s-i :  (lnxB)^i * (lnQ2)^j
//   Block B:  for k=1..dR, for m=0..dZQ      :  (lnQ2)^m * z^k
//   Block C:  for k=1..dR, for n=1..dZX      :  (lnxB)^n * z^k
void correctionTools::buildMcBasis(double xB, double Q2, double z,
                                   int dP, int dR, int dZQ, int dZX,
                                   std::vector<double>& f){
	const int M = mcBasisSize(dP, dR, dZQ, dZX);
	f.assign(M, 0.0);

	// Guard against pathological inputs to log()
	const double lnX  = (xB > 0) ? std::log(xB) : 0.0;
	const double lnQ2 = (Q2 > 0) ? std::log(Q2) : 0.0;

	// Precompute powers
	std::vector<double> powLnX (dP + dZX + 1, 1.0);
	std::vector<double> powLnQ2(dP + dZQ + 1, 1.0);
	for( int i = 1; i < (int)powLnX.size();  i++ ) powLnX[i]  = powLnX[i-1]  * lnX;
	for( int i = 1; i < (int)powLnQ2.size(); i++ ) powLnQ2[i] = powLnQ2[i-1] * lnQ2;
	std::vector<double> powZ(dR + 1, 1.0);
	for( int k = 1; k <= dR; k++ ) powZ[k] = powZ[k-1] * z;

	int idx = 0;

	// ── Block A: polynomial in (lnxB, lnQ2), total degree <= dP ──
	for( int s = 0; s <= dP; s++ ){
		for( int i = 0; i <= s; i++ ){
			int j = s - i;
			f[idx++] = powLnX[i] * powLnQ2[j];
		}
	}

	// ── Block B: z^k * (lnQ2)^m, k=1..dR, m=0..dZQ ──
	for( int k = 1; k <= dR; k++ ){
		for( int m = 0; m <= dZQ; m++ ){
			f[idx++] = powLnQ2[m] * powZ[k];
		}
	}

	// ── Block C: z^k * (lnxB)^n, k=1..dR, n=1..dZX ──
	for( int k = 1; k <= dR; k++ ){
		for( int n = 1; n <= dZX; n++ ){
			f[idx++] = powLnX[n] * powZ[k];
		}
	}

	// idx must equal M; this is asserted-by-construction.
	(void)M; // suppress unused warning in non-debug builds
}

void correctionTools::loadFits(){
	
	TFile * weightFile = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	TFile * pi2kFile = new TFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	TFile * k2piFile = new TFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	
	for( int x = 0; x < bins_xB; x++ ){
		for( int q = 0; q < bins_Q2; q++ ){
			weightFit[0][x][q] = (TF1*)weightFile->Get( Form("fitPip_%i_%i", q, x) );
			weightFit[1][x][q] = (TF1*)weightFile->Get( Form("fitPim_%i_%i", q, x) );

			for( int p = 0; p < 4; p++ ){
				k2piFit[0][x][q][p] = (TF1*)k2piFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				k2piFit[1][x][q][p] = (TF1*)k2piFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				
				pi2kFit[0][x][q][p] = (TF1*)pi2kFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
				pi2kFit[1][x][q][p] = (TF1*)pi2kFile->Get( Form("fitPip_%i_%i_%i", p, q, x) );
			}
		}
	}
	weightFile->Close();
	pi2kFile->Close();
	k2piFile->Close();
}

void correctionTools::loadParameters(){
	printFilePaths();
	loadWeightParameters();
	std::cout<<"Loaded Weight Parameters\n";
	loadkTopiParameters();
	std::cout<<"Loaded k-to-pi Parameters\n";
	loadpiTokParameters();
	//std::cout<<"Loaded pi-to-k Parameters\n";
}	

void correctionTools::loadWeightParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	std::cout<<"Opened weight file\n";
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 3; j++ ){
			
			TF1 * fp = (TF1*)weightFile.Get( Form("fitPipQ_f_%i_%i", i, j));
			weight_Parameters[0][i][j][0] = fp->GetParameter(0);
			weight_Parameters[0][i][j][1] = fp->GetParameter(1);
			weight_Parameters[0][i][j][2] = fp->GetParameter(2);
			weight_Parameters[0][i][j][3] = fp->GetParameter(3);
			
			TF1 * fm = (TF1*)weightFile.Get(Form("fitPimQ_f_%i_%i", i, j));
			weight_Parameters[1][i][j][0] = fm->GetParameter(0);
			weight_Parameters[1][i][j][1] = fm->GetParameter(1);
			weight_Parameters[1][i][j][2] = fm->GetParameter(2);
			weight_Parameters[1][i][j][3] = fm->GetParameter(3);
		}
	}

}

void correctionTools::loadkTopiParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);
	std::cout<<"Opened p2k file\n";
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 4; j++ ){
			for( int k = 2; k < 4; k++ ){
				TF1 * fp = (TF1*)weightFile.Get(Form("fitPipX_%i_%i_%i", k, i, j));
				k_To_pi_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				k_To_pi_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				k_To_pi_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				k_To_pi_Parameters[0][k][i][j][3] = fp->GetParameter(3);
				//k_To_pi_Parameters[0][k][i][j][4] = fp->GetParameter(4);
				
				TF1 * fm = (TF1*)weightFile.Get(Form("fitPimX_%i_%i_%i", k, i, j));
				k_To_pi_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				k_To_pi_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				k_To_pi_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				k_To_pi_Parameters[1][k][i][j][3] = fm->GetParameter(3);
				//k_To_pi_Parameters[1][k][i][j][4] = fm->GetParameter(4);
				
				//std::cout<<"Filled Parameters : i = "<<i<<" , j = "<<j<<std::endl;
			}
		}
	}
	
}

void correctionTools::loadpiTokParameters(){
	TFile weightFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	for( int i = 0; i < 4; i++ ){//momentum bins
		for( int j = 0; j < 4; j++ ){//loop through fits for each a(x, q)
			for( int k = 2; k < 4; k++ ){//loop through fits for each b(x)
				TF1 * fp = (TF1*)weightFile.Get(Form("fitPipX_%i_%i_%i", k, i, j));
				pi_To_k_Parameters[0][k][i][j][0] = fp->GetParameter(0);
				pi_To_k_Parameters[0][k][i][j][1] = fp->GetParameter(1);
				pi_To_k_Parameters[0][k][i][j][2] = fp->GetParameter(2);
				pi_To_k_Parameters[0][k][i][j][3] = fp->GetParameter(3);
				//pi_To_k_Parameters[0][k][i][j][4] = fp->GetParameter(4);
			
				TF1 * fm = (TF1*)weightFile.Get(Form("fitPimX_%i_%i_%i", k, i, j));
				pi_To_k_Parameters[1][k][i][j][0] = fm->GetParameter(0);
				pi_To_k_Parameters[1][k][i][j][1] = fm->GetParameter(1);
				pi_To_k_Parameters[1][k][i][j][2] = fm->GetParameter(2);
				pi_To_k_Parameters[1][k][i][j][3] = fm->GetParameter(3);
				//pi_To_k_Parameters[1][k][i][j][4] = fm->GetParameter(4);
			
			}
		}
	}
}

void correctionTools::setKinematics( double x, double q, double z, double p ){ 
	kin[0] = x;
	kin[1] = q;
	kin[2] = z;
	kin[3] = p;
}


double correctionTools::getCorrectionFactor( int type, int charge ){
	if( charge != 0 && charge != 1 ){
		std::cout<<"Invalid charge\n";
		return 0;
	}

	double x = kin[0];
	double q = kin[1];
	double z = kin[2];
	double p = kin[3];

	// ── Continuous fit overrides (loaded by loadContinuousCorrections) ────────
	// When mode==4 and per-bin 4D fits are loaded, prefer them over the integrated 3D fit.
	// bin_4d==-1 means the caller did not select a bin (out-of-range event); fall back to 3D.
	if( type == 2 && mode == 4 && nBins_4d > 0 && bin_4d >= 0 ){
		int b = (bin_4d < nBins_4d) ? bin_4d : nBins_4d - 1;
		if( mcCorrectionFit4D[charge][b] )
			return mcCorrectionFit4D[charge][b]->Eval(x, q, z);
	}
	if( type == 2 && mcCorrectionFit[charge] )
		return mcCorrectionFit[charge]->Eval(x, q, z);
	if( type == 2 && nBins_4d > 0 ){
		int b = (bin_4d >= 0 && bin_4d < nBins_4d) ? bin_4d : 0;
		if( mcCorrectionFit4D[charge][b] )
			return mcCorrectionFit4D[charge][b]->Eval(x, q, z);
		return 1.0;
	}
	if( type == 3 && pi2kRatioFit[charge] )
		return (p >= misid_min) ? pi2kRatioFit[charge]->Eval(p) : 1.0;
	if( type == 4 && k2piRatioFit[charge] )
		return (p >= misid_min) ? k2piRatioFit[charge]->Eval(p) : 0.0;
	if( type == 3 || type == 4 ){
		if( !pi_to_k_Correction[charge] || !k_to_pi_Correction[charge] ) return (type == 3) ? 1.0 : 0.0;
		int p_bin = pi_to_k_Correction[0]->FindBin(p);
		if( type == 3 ) return (p >= misid_min) ? pi_to_k_Correction[charge]->GetBinContent(p_bin) : 1.0;
		if( type == 4 ) return (p >= misid_min) ? k_to_pi_Correction[charge]->GetBinContent(p_bin) : 0.0;
	}
	// ── TH3 binned path (only reached when loadHistograms was called) ─────────
	if( !binMigration[charge] || !accCorrection[charge] ) return 1.0;

	int x_bin_w = binMigration[charge]->GetXaxis()->FindBin(x);
	int y_bin_w = binMigration[charge]->GetYaxis()->FindBin(q);
	int z_bin_w = binMigration[charge]->GetZaxis()->FindBin(z);
	double binWeight = binMigration[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
	double accWeight = accCorrection[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w);
	double mcWeight  = mcCorrection[charge]  ? mcCorrection[charge]->GetBinContent( x_bin_w, y_bin_w, z_bin_w) : 1.0;

	if( type == 0 ) return binWeight;
	if( type == 1 ) return accWeight;
	if( type == 2 ) return mcWeight;
	if( type > 4 || type < 0 ){ return 0; }	

	/*
	//Check if in correction phase space       
	double badBinVals[2] = {1, 0}; //Don't correct "bad" pi2k bins, kill "bad" k2pi bins

	if( kWeight <= 0 ){return badBinVals[type - 3];}  //|| kWeight > 1 || this_bin_p < 2
	else{
		double kWeightTemp;
		int tempCharge;
		if( charge == 0 ){ tempCharge = 1; }
		if( charge == 1 ){ tempCharge = 0; }

		if(type == 3){
			kWeightTemp = pi_to_k_Correction[tempCharge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
		}
		if( type == 4 ){
			kWeightTemp = k_to_pi_Correction[tempCharge][this_bin_p]->GetBinContent( x_bin_k, y_bin_k, z_bin_k );
		}

		if( kWeightTemp <= 0 || kWeightTemp > 1 ){
			return badBinVals[type-3];
		}
		return kWeight;
	}
	*/
	return 1;

	

}		

double correctionTools::getCorrectionError( int type, int charge ){
	if( charge != 0 && charge != 1 ){
		std::cout<<"Invalid charge\n";
		return 0;
	}

	double x = kin[0];
	double q = kin[1];
	double z = kin[2];
	double p = kin[3];

	// ── Continuous fit overrides (loaded by loadContinuousCorrections) ────────
	if( type == 2 && mcCorrectionFitSigma[charge] )
		return mcCorrectionFitSigma[charge]->Eval(x, q, z);
	if( type == 2 && nBins_4d > 0 ){
		int b = (bin_4d >= 0 && bin_4d < nBins_4d) ? bin_4d : 0;
		if( mcCorrectionFit4DSigma[charge][b] )
			return mcCorrectionFit4DSigma[charge][b]->Eval(x, q, z);
		return 0.0;
	}
	if( type == 3 && pi2kRatioFitSigma[charge] )
		return (p >= misid_min) ? pi2kRatioFitSigma[charge]->Eval(p) : 0.0;
	if( type == 4 && k2piRatioFitSigma[charge] )
		return (p >= misid_min) ? k2piRatioFitSigma[charge]->Eval(p) : 0.0;
	// ── TH3 binned path (only reached when loadHistograms was called) ─────────
	if( !binMigration[charge] ) return 0.0;

	int x_bin_w = binMigration[charge]->GetXaxis()->FindBin(x);
	int y_bin_w = binMigration[charge]->GetYaxis()->FindBin(q);
	int z_bin_w = binMigration[charge]->GetZaxis()->FindBin(z);

	double binErr  = binMigration[charge]->GetBinError( x_bin_w, y_bin_w, z_bin_w);
	double accErr  = accCorrection[charge] ? accCorrection[charge]->GetBinError( x_bin_w, y_bin_w, z_bin_w) : 0.0;
	double mcErr   = mcCorrection[charge]  ? mcCorrection[charge]->GetBinError(  x_bin_w, y_bin_w, z_bin_w) : 0.0;
	double pi2kErr = 0.0, k2piErr = 0.0;
	if( pi_to_k_Correction[charge] && k_to_pi_Correction[charge] ){
		int p_bin = pi_to_k_Correction[charge]->GetXaxis()->FindBin(p);
		pi2kErr = (p >= misid_min) ? pi_to_k_Correction[charge]->GetBinError(p_bin) : 0.0;
		k2piErr = (p >= misid_min) ? k_to_pi_Correction[charge]->GetBinError(p_bin) : 0.0;
	}

	if( type == 0 ) return binErr;
	if( type == 1 ) return accErr;
	if( type == 2 ) return mcErr;
	if( type == 3 ) return pi2kErr;
	if( type == 4 ) return k2piErr;

	return 0;
}		
