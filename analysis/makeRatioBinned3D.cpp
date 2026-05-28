#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <regex>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom3.h"
#include "TParameter.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"
#include "analyzer.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST


using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;
using std::isnan;
using namespace cutVals;
using namespace constants;
// -------------------------------
// Type aliases for readability
// -------------------------------
using VecE      = std::vector<double>;   // size N_E
using VecP      = std::vector<VecE>;     // size bins_p
using VecZ      = std::vector<VecP>;     // size bins_Z
using VecXB     = std::vector<VecZ>;     // size bins_xB
using VecQ2     = std::vector<VecXB>;    // size bins_Q2
using VecVar    = std::vector<VecQ2>;    // size nVar
using VecCharge = std::vector<VecVar>;   // size nCharge
using VecSample = std::vector<VecCharge>;// size nSamples

// For weights: no p-dimension
using VecZ_w    = std::vector<VecE>;     
using VecXB_w   = std::vector<VecZ_w>;
using VecQ2_w   = std::vector<VecXB_w>;
using VecVar_w  = std::vector<VecQ2_w>;
using VecCharge_w = std::vector<VecVar_w>;
using VecSample_w = std::vector<VecCharge_w>;

//config... edited by input
int matchType = 0;
int applyCorr = 0;
int binnedOut = 1;
int map = 0;
int n4d_corr_bins = 0;  // number of 4D correction bins found in file; 0 = use 3D corrections
int target = 0;         // 0 = RGB/deuterium, 1 = RGA/proton
TString var_name = "null";
int bins_var = 1;
double var_min = 0;
double var_max = 1;
TString matching_file = "matchCut2D_10.root";
int overflowBinning = 0;  // bitmap: 1=Q2, 2=xB, 4=Z, 8=var  (e.g. 7=all kinematic, 15=all)
int verbose = 0;          // 0 = concise output; 1 = per-bin diagnostic breakdown of Var(sigma+/-)

void zeroSuppress( TH1F * h);
double getVarVal(TString var,  electron e, pion pi );
void setBin(
    TH1F* h,
    int z_bin,
    const VecZ& events,      // [Z][p][E]
    const VecZ& errors,      // [Z][p][E]
    double weights,
    correctionTools& corr,
	int corrType,
    int chargeIdx,
    int E_bin,
    TString corrName
);

bool updateCorrectionsForBeam(
    double eBeam,
    double& beam_energy,
	TString& corrName,
    correctionTools& corrector
);
bool updateCorrectionsForBeam(
	int E,
    TString& corrName,
    correctionTools& corrector
);
void getNewRhoName( double E, TString& rho_norm_name) ;


int main( int argc, char** argv){

	if( argc < 13 ){
		cerr << "Usage: ./makeRatioBinned3D [beam energy] [output file]\n";
		cerr << "       [match type] [apply corr] [map] [acceptance file]\n"
	     << "       apply corr: bitwise — bits[1:0]: 0=none,1=bin migration,2=acceptance,3=full MC; bit2(4)=kaon k-factor; bit3(8)=rho/VM\n";
		cerr << "       [var name] [bins var] [var min] [var max] [binned out] [correction name]\n";
		cerr << "       (opt: Mx2 cut) (opt: rho mode) (opt: target) (opt: matching file)\n";
		cerr << "       (opt: overflow binning) (opt: verbose 0/1)\n";
		return -1;
	}

	double inBeam = atof(argv[1]);

	TString out_name = argv[2];
	matchType = atoi(argv[3]);
	applyCorr = atoi(argv[4]);
	map = atoi(argv[5]);
	TString acc_name = argv[6];

	var_name = argv[7];
	bins_var = atoi(argv[8]);
	var_min = atof(argv[9]);
	var_max = atof(argv[10]);
	if( atoi(argv[11]) > 1 ){
		binnedOut = bins_var;
	}
	TString correction_name = argv[12];
	double Mx2_cut = 1.25;
	double rho_mode = 1;
	if( argc > 13 ) Mx2_cut = atof(argv[13]);
	if( argc > 14 ) rho_mode = atoi(argv[14]);
	if( argc > 15 ) target = atoi(argv[15]);
	if( argc > 16) matching_file = argv[16];
	if( argc > 17) overflowBinning = atoi(argv[17]);
	if( argc > 18) verbose = atoi(argv[18]);

	cout << "Target: " << target << " (" << (target == 1 ? "RGA/proton" : "RGB/deuterium") << ")\n";
	if (overflowBinning)
		cout << "Overflow binning bitmap=" << overflowBinning
			 << " (" << ((overflowBinning & 1) ? "Q2 " : "")
			 << ((overflowBinning & 2) ? "xB " : "")
			 << ((overflowBinning & 4) ? "Z "  : "")
			 << ((overflowBinning & 8) ? "var" : "") << ")\n";

	// Probe correction file for 4D binning metadata before any array allocation
	if( bins_var > 1 ){
		TString probePath = (TString)_DATA + "/correctionFiles/" + correction_name;
		TFile * pf = TFile::Open(probePath);
		if( pf && !pf->IsZombie() ){
			int count = 0;
			while( pf->Get(Form("hMcCorrP_fit_%d", count)) ) count++;
			n4d_corr_bins = count;
			if( count > 0 ){
				auto* pMin = (TParameter<double>*)pf->Get("var_min");
				auto* pMax = (TParameter<double>*)pf->Get("var_max");
				if( pMin ) var_min = pMin->GetVal();
				if( pMax ) var_max = pMax->GetVal();
				cout << "Correction file: " << count << " 4D bins, var=[" << var_min << ", " << var_max << "]\n";
			} else {
				cout << "No 4D histograms in correction file; using 3D corrections (xB, Q2, z)\n";
			}
			delete pf;
		} else {
			cerr << "Warning: could not open correction file for probing: " << probePath << "\n";
		}
	}

	cout << "Beam energy: " << inBeam << "  Output: " << out_name << "\n";
	cout << "Match type: " << matchType << "  Apply corr: " << applyCorr << "  Map: " << map << "\n";
	cout << "Var: " << var_name << " [" << var_min << ", " << var_max << "] (" << bins_var << " bins)\n";
	cout << "Corrections: " << correction_name << "  Mx2 cut: " << Mx2_cut << "\n";


	TFile * outFile = new TFile(out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TChain * dChain = new TChain( "ePi" ); 
	TChain * kChain = new TChain( "ePi" );
	TChain * rChain = new TChain( "ePi" );

	//Add files
	TString base = "../trees/"; //"/volatile/clas12/users/jphelan/SIDIS/data/";
	constexpr int N_E = 3;
	

	if( target == 0 ){
		if( inBeam == 0 || inBeam == 10.2){
			dChain->Add( base + "final_skims/10.2/final_skim.root");
			if( applyCorr & 4 ){
				kChain->Add( base + "final_skims/kaons_10.2/final_skim.root");
			}
			if( applyCorr & 8 ){
				rChain->Add( base + "final_skims/rho_skims/rotated_10.2_final.root");
			}
		}
		if( inBeam == 0 || inBeam == 10.4){
			dChain->Add( base + "final_skims/10.4/final_skim.root");
			if( applyCorr & 4 ){
				kChain->Add( base + "final_skims/kaons_10.4/final_skim.root");
			}
			if( applyCorr & 8 ){
				rChain->Add( base + "final_skims/rho_skims/rotated_10.4_final.root");
			}
		}
		if( inBeam == 0 || inBeam == 10.6){
			dChain->Add( base + "final_skims/10.6/final_skim.root");
			if( applyCorr & 4 ){
				kChain->Add( base + "final_skims/kaons_10.6/final_skim.root");
			}
			if( applyCorr & 8 ){
				rChain->Add( base + "final_skims/rho_skims/rotated_10.6_final.root");
			}
		}
	}
	if( target == 1 ){
		dChain->Add( base + "final_skims/10.2/final_skim_rga.root");
		if( applyCorr & 4 ){
			kChain->Add( base + "final_skims/kaons_10.2/final_skim_rga.root");
		}
		if( applyCorr & 8 ){
			rChain->Add( "../trees/rga/rotated_10.2_final.root");
		}
	}



	cout<<"Chains loaded\n";



// -------------------------------
// Dimensions
// -------------------------------
constexpr int nSamples = 4;    // sample
constexpr int nCharge  = 2;    // charge
int nVar = bins_var;           // number of variables (runtime)


// -------------------------------
// Events in bin: [sample][charge][var][Q2][xB][Z][p][E]
// -------------------------------
VecSample events_in_bin(
    nSamples, VecCharge(
        nCharge, VecVar(
            nVar, VecQ2(
                bins_Q2, VecXB(
                    bins_xB, VecZ(
                        bins_Z, VecP(
                            bins_p, VecE(N_E, 0.0)
                        )
                    )
                )
            )
        )
    )
);

// -------------------------------
// Errors in bin: [sample][charge][var][Q2][xB][Z][p][E]
// -------------------------------
VecSample errors_in_bin(
    nSamples, VecCharge(
        nCharge, VecVar(
            nVar, VecQ2(
                bins_Q2, VecXB(
                    bins_xB, VecZ(
                        bins_Z, VecP(
                            bins_p, VecE(N_E, 0.0)
                        )
                    )
                )
            )
        )
    )
);

VecSample_w weights_in_bin(
    nSamples, VecCharge_w(
        nCharge, VecVar_w(
            nVar, VecQ2_w(
                bins_Q2, VecXB_w(
                    bins_xB, VecZ_w(
                        bins_Z, VecE(N_E, 0.0)
                    )
                )
            )
        )
    )
);

// -------------------------------
// VM weight uncertainty sum: [charge][var][Q2][xB][Z][E]
//
// Stores Var(W^h) = sum_i (delta_alpha_i^h)^2  (formula (2)).
// Filled only by the rho/VM loop.  Since the VM subtraction is
// integrated over p (alpha does not depend on a p-bin), no p
// dimension is needed here.  This accumulator holds *only* the
// weight-uncertainty piece — the Poisson counting variance of the
// alpha_i themselves is not included, matching the convention that
// Var(W^h) is the systematic uncertainty on the VM weights alone.
// -------------------------------
VecCharge_w vm_weight_var_in_bin(
    nCharge, VecVar_w(
        nVar, VecQ2_w(
            bins_Q2, VecXB_w(
                bins_xB, VecZ_w(
                    bins_Z, VecE(N_E, 0.0)
                )
            )
        )
    )
);

// -------------------------------
// RICH-decomposed counts: [sample][charge][var][Q2][xB][Z][p][E]
// Tracks events in RICH acceptance for pion (sample 0) and kaon (sample 1) samples
// Used for full Poisson uncertainty propagation of k-factors
// -------------------------------
VecSample rich_in_bin(
    nSamples, VecCharge(
        nCharge, VecVar(
            nVar, VecQ2(
                bins_Q2, VecXB(
                    bins_xB, VecZ(
                        bins_Z, VecP(
                            bins_p, VecE(N_E, 0.0)
                        )
                    )
                )
            )
        )
    )
);

// -------------------------------
// Event-weighted kinematic sums: [var][Q2][xB][Z][E]
// Used to compute mean kinematics per overflow bin for correction lookup
// -------------------------------
VecVar_w sum_Q2_kin(
    nVar, VecQ2_w(
        bins_Q2, VecXB_w(
            bins_xB, VecZ_w(
                bins_Z, VecE(N_E, 0.0)
            )
        )
    )
);
VecVar_w sum_xB_kin(
    nVar, VecQ2_w(
        bins_Q2, VecXB_w(
            bins_xB, VecZ_w(
                bins_Z, VecE(N_E, 0.0)
            )
        )
    )
);
VecVar_w sum_z_kin(
    nVar, VecQ2_w(
        bins_Q2, VecXB_w(
            bins_xB, VecZ_w(
                bins_Z, VecE(N_E, 0.0)
            )
        )
    )
);
VecVar_w count_kin(
    nVar, VecQ2_w(
        bins_Q2, VecXB_w(
            bins_xB, VecZ_w(
                bins_Z, VecE(N_E, 0.0)
            )
        )
    )
);


	TH1F * hZ[bins_var][bins_Q2][bins_xB][2];
	TH1F * hZ_k[bins_var][bins_Q2][bins_xB][2];
	TH1F * hZ_r[bins_var][bins_Q2][bins_xB][2];

	TString charge_str[2] = {"", "_Pim"};

	cout<<"Make histograms\n";
	for (int var = 0; var < bins_var; var++) {
		for (int i = 0; i < bins_Q2; i++) {
			for (int j = 0; j < bins_xB; j++) {
				for (int k = 0; k < 2; k++) {
	
					TString suffix;
	
					if (bins_var > 1)
						suffix = Form("_%d_%d_%d", i+1, j+1, k+1);
					else
						suffix = Form("_%d_%d", i+1, j+1);
	
					// Base piece reused everywhere
					TString base = Form("_%d%s%s",
										var,
										charge_str[k].Data(),
										suffix.Data());
	
					// ---- hRatio ----
					hZ[var][i][j][k] = new TH1F(
						"hRatio" + base,
						Form("hRatio %s%s",
							 charge_str[k].Data(),
							 suffix.Data()),
						bins_Z, 0.3, 1.0
					);
	
					hZ[var][i][j][k]->Sumw2();

					// ---- hRatio_k ----
					hZ_k[var][i][j][k] = new TH1F(
						"hRatio_k" + base,
						Form("hRatio_k %s%s",
							 charge_str[k].Data(),
							 suffix.Data()),
						bins_Z, 0.3, 1.0
					);
	
					hZ_k[var][i][j][k]->Sumw2();

					// ---- hRatio_r ----
					hZ_r[var][i][j][k] = new TH1F(
						"hRatio_r" + base,
						Form("hRatio_r %s%s",
							 charge_str[k].Data(),
							 suffix.Data()),
						bins_Z, 0.3, 1.0
					);
	
					hZ_r[var][i][j][k]->Sumw2();

				
		
				}
			}
		}
	}

	correctionTools corrector(n4d_corr_bins > 0 ? 4 : 2);
	//corrector.setWeightName( correction_name );
	//corrector.setWeightName4D( correction_name );
	if( n4d_corr_bins > 0 ){
		corrector.setN4dBins(n4d_corr_bins);
	}
	//corrector.setK2piName( target == 1 ? "corrections_misid_rga.root" : "corrections_misid_AN.root");
	//corrector.setPi2kName( target == 1 ? "corrections_misid_rga.root" : "corrections_misid_AN.root");

	corrector.setWeightFitName( correction_name );
	corrector.setMisidFitName(  "corrections_misid_fit.root" );

	//corrector.loadHistograms();	
	corrector.loadContinuousCorrections();
	corrector.printFilePaths();

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions(matching_file);
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/"+acc_name);//%.1f.root", energy));

	////////////////////////////
	///////// Pions ////////////
	////////////////////////////

	cout << "Begin pion analysis\n";
	double beam_energy = 10.2; //current energy of file

	//TTreeReader reader_rec("ePi", inFile);
	TTreeReader reader_rec(dChain);
	TTreeReaderValue<double> eBeam( reader_rec, "Ebeam" );
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");

	TTreeReaderArray<bool> isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion3D(reader_rec, "isGoodPion_3d");

	// Lambda: compute all bin indices for a given electron/pion/momentum
	struct Bins { int E, Q2, xB, Z, var, p; };
	auto calcBins = [&](const electron& el, const pion& pi_ref, double p_mom) -> Bins {
		int bE   = (int)((beam_energy - 10.2) / .2);
		int bQ2  = (int)(((el.getQ2()  - Q2_min) / (Q2_max - Q2_min)) * bins_Q2);
		int bxB  = (int)(((el.getXb()  - xB_min) / (xB_max - xB_min)) * bins_xB);
		int bZ   = (int)(((pi_ref.getZ() - .3)    / (1. - .3))         * bins_Z);
		int bVar = (int)(((getVarVal(var_name, el, pi_ref) - var_min) / (var_max - var_min)) * bins_var);
		int bP   = -1;
		for (int bin = 0; bin < bins_p; bin++)
			if (p_mom > p_bin_edges[bin] && p_mom < p_bin_edges[bin+1]) { bP = bin; break; }

		// Safety: reject under-range events (negative bin index = invalid).
		// For variables with overflow bit set, events above max are valid
		// (they'll be filled into all bins >= their native bin).
		// For variables without overflow, reject events outside [min, max).
		if (bQ2  < 0 || (!(overflowBinning & 1) && bQ2  >= bins_Q2))  bQ2  = -1;
		if (bxB  < 0 || (!(overflowBinning & 2) && bxB  >= bins_xB))  bxB  = -1;
		if (bZ   < 0 || (!(overflowBinning & 4) && bZ   >= bins_Z))    bZ   = -1;
		if (bVar < 0 || (!(overflowBinning & 8) && bVar >= bins_var))  bVar = -1;
		if (bE   < 0 || bE  >= N_E)     bE  = -1;
		// Clamp overflow events to the last physical bin (their "minimum" bin)
		if ((overflowBinning & 1) && bQ2  >= bins_Q2)  bQ2  = bins_Q2  - 1;
		if ((overflowBinning & 2) && bxB  >= bins_xB)  bxB  = bins_xB  - 1;
		if ((overflowBinning & 4) && bZ   >= bins_Z)    bZ   = bins_Z   - 1;
		if ((overflowBinning & 8) && bVar >= bins_var)  bVar = bins_var - 1;

		return {bE, bQ2, bxB, bZ, bVar, bP};
	};

	// Lambda: apply acceptance matching based on matchType
	auto applyMatching = [&](const pion& pi_ref, double p_mom) -> bool {
		int chargeIdx = (int)( pi_ref.getCharge() < 1 )+1;
		if (matchType == 1 || matchType == 2)
			return anal.applyAcceptanceMatching(pi_ref, 2);
		int chargeIdx_1 = ( matchType == 3 ) ? 1 : chargeIdx;
		int chargeIdx_2 = ( matchType == 3 ) ? 2 : chargeIdx;

		if ( matchType == 3 || matchType == 4 )
			return  anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 1, chargeIdx_1) >= 0 &&
			       anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 2, chargeIdx_2) >= 0 && anal.applyAcceptanceMatching(pi_ref, 2);
		return true;
	};

	int event_total = reader_rec.GetEntries();
	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();
		if(event_count%100000 == 0){
			cout << "Pions: " << event_count << " / " << event_total << "\n";
		}

		updateCorrectionsForBeam(*eBeam, beam_energy, correction_name, corrector);

		TVector3 e_mom = e->get3Momentum();
		if(matchType<=1 && map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			int pid = abs(pi[i].getPID_eb());
			if( pid > 321 ) continue;
			int pid_bin = (pid==211) ? 0 : 1;
			double p_pi = pi[i].get3Momentum().Mag();

			if(matchType==0 && !isGoodPion_no_acc[i] ) continue;
			else if(matchType==1 && !isGoodPion_no_acc[i] ) continue;
			else if(matchType==2 && !isGoodPion_no_acc[i]) continue;
			else if(matchType==3 && (!isGoodPion_no_acc[i])) continue;
			else if(matchType==4 && (!isGoodPion_no_acc[i])) continue;

			if(matchType<=1 && map && anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

			if( !applyMatching(pi[i], p_pi) ) continue;

			auto b = calcBins(*e, pi[i], p_pi);
			if (b.E < 0 || b.Q2 < 0 || b.xB < 0 || b.Z < 0 || b.var < 0) continue;

			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );
			double weight = 1;
			
			if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
			if( (applyCorr & 3) == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( (applyCorr & 3) == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( (applyCorr & 3) == 3 ) weight *= corrector.getCorrectionFactor(2, chargeIdx);
			if( (applyCorr & 4) && b.p >= 0) weight *= (pid==211) ? corrector.getCorrectionFactor(3, chargeIdx) : corrector.getCorrectionFactor(4, chargeIdx) ;
			
			// Overflow binning: every event fills ALL bins for flagged variables.
			// Each overflow bin is identical (integrated over that variable).
			// Bounded mode: fill only the native bin.
			int Q2_lo = (overflowBinning & 1) ? 0     : b.Q2;
			int Q2_hi = b.Q2 + 1;
			int xB_lo = (overflowBinning & 2) ? 0     : b.xB;
			int xB_hi = b.xB + 1;
			int Z_lo  = (overflowBinning & 4) ? 0     : b.Z;
			int Z_hi  = b.Z  + 1;
			for (int iq = Q2_lo; iq < Q2_hi; iq++) {
				for (int ix = xB_lo; ix < xB_hi; ix++) {
					for (int iz = Z_lo; iz < Z_hi; iz++) {
						if( !(applyCorr & 4) && pid!=211 ) continue;
						events_in_bin[pid_bin][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;
						weights_in_bin[pid_bin][chargeIdx][b.var][iq][ix][iz][b.E] += weight;
						errors_in_bin[pid_bin][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;

						// Track RICH acceptance for full Poisson uncertainty on k-factors
						if (b.p >= 0 && pi[i].getBeta_rich() > 0) {
							rich_in_bin[pid_bin][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;
						}

						// Accumulate event kinematics for mean correction lookup
						sum_Q2_kin[b.var][iq][ix][iz][b.E] += e->getQ2();
						sum_xB_kin[b.var][iq][ix][iz][b.E] += e->getXb();
						sum_z_kin [b.var][iq][ix][iz][b.E] += pi[i].getZ();
						count_kin [b.var][iq][ix][iz][b.E] += 1.0;
					}
				}
			}

		}
	}
/*
	///////////////////////////////////////
	////////////// Kaons /////////////////
	///////////////////////////////////////

	//TTreeReader reader_k("ePi", kFile);
	TTreeReader reader_k(kChain);
	TTreeReaderValue<double> eBeam_k( reader_k, "Ebeam" );
	TTreeReaderValue<electron> e_k(reader_k, "e");
	TTreeReaderArray<pion> k(reader_k, "pi");
	TTreeReaderArray<bool> isGoodKaon_no_acc(reader_k, "isGoodPion_no_acc");
	TTreeReaderArray<bool> isGoodKaon(reader_k, "isGoodPion");
	TTreeReaderArray<bool> isGoodKaon3D(reader_k, "isGoodPion_3d");


	event_total = reader_k.GetEntries();

	if( applyCorr & 4 ){
		cout << "Begin kaon analysis\n";
		while (reader_k.Next()) {
			int event_count = reader_k.GetCurrentEntry();
			if(event_count%100000 == 0){
				cout << "Kaons: " << event_count << " / " << event_total << "\n";
			}
			
			updateCorrectionsForBeam(*eBeam_k, beam_energy, correction_name, corrector);

			TVector3 e_mom = e_k->get3Momentum();
			if(matchType==0 && map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
	
			for( int i = 0; i < (int) ( k.end() - k.begin() ); i++ ){
				
				int chargeIdx = (int)( k[i].getCharge() < 1 );
				double p_pi = k[i].get3Momentum().Mag();
				if( matchType == 0 && !isGoodKaon_no_acc[i] )continue;
				else if( matchType == 2 && !isGoodKaon[i] )continue;
				else if( matchType == 3 && !isGoodKaon3D[i] )continue;

				if(matchType==0 && map && anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

				//if( !applyMatching(k[i], p_pi) ) continue;

				auto b = calcBins(*e_k, k[i], p_pi);
				if (b.E < 0 || b.Q2 < 0 || b.xB < 0 || b.Z < 0 || b.var < 0) continue;
				corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );
				if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
				double weight = 1;
				if( (applyCorr & 3) == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
				if( (applyCorr & 3) == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
				if( (applyCorr & 3) == 3 ) weight *= corrector.getCorrectionFactor(2, chargeIdx);
				if( applyCorr & 4 )        weight *= corrector.getCorrectionFactor(4, chargeIdx);
				int Q2_lo = (overflowBinning & 1) ? 0     : b.Q2;
				int Q2_hi = b.Q2 + 1;
				int xB_lo = (overflowBinning & 2) ? 0     : b.xB;
				int xB_hi = b.xB + 1;
				int Z_lo  = (overflowBinning & 4) ? 0     : b.Z;
				int Z_hi  = b.Z  + 1;
				for (int iq = Q2_lo; iq < Q2_hi; iq++) {
					for (int ix = xB_lo; ix < xB_hi; ix++) {
						for (int iz = Z_lo; iz < Z_hi; iz++) {
							events_in_bin[1][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;
							weights_in_bin[1][chargeIdx][b.var][iq][ix][iz][b.E] += weight;
							errors_in_bin[1][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;

							// Track RICH acceptance for full Poisson uncertainty on k-factors
							if (b.p >= 0 && k[i].getBeta_rich() > 0) {
								rich_in_bin[1][chargeIdx][b.var][iq][ix][iz][b.p][b.E]++;
							}
						}
					}
				}

			}
		}
	}
*/
	///////////////////////////////////////
	////////////// Rhos /////////////////
	///////////////////////////////////////
	
	//TTreeReader reader_r("ePi", rFile);
	TTreeReader reader_r(rChain);

	TTreeReaderValue<electron> e_r(reader_r, "e");
	TTreeReaderArray<pion> r(reader_r, "pi");
	TTreeReaderArray<double> rhoWeight( reader_r, "rhoWeight");
	TTreeReaderValue<double> Mx_2pi( reader_r, "Mx_2pi");
	TTreeReaderArray<double> rhoError( reader_r, "rhoErr");
	//TTreeReaderArray<double> rhoError_sym( reader_r, "rhoErr_sym");
	TTreeReaderArray<bool> isGoodRho(reader_r, "isGoodPion");
	TTreeReaderValue<TLorentzVector> beam(reader_r, "beam");

	event_total = reader_r.GetEntries();

	if( applyCorr & 8 ){
		cout << "Begin rho analysis\n";
		while (reader_r.Next()) {
			int event_count = reader_r.GetCurrentEntry();

			if(event_count%100000 == 0){
				cout << "Rhos:  " << event_count << " / " << event_total << "\n";
			}
			updateCorrectionsForBeam(beam->E(), beam_energy, correction_name, corrector);
			
			
			TVector3 e_mom = e_r->get3Momentum();

			int i = -1; 
			
			// If our analysis pion is good, continue with that pion.  If our analysis pion is not good, 
			if      ( isGoodRho[rho_mode] ) i = rho_mode;
			else if ( rho_mode==1 && !isGoodRho[1] && r[1].get3Momentum().Theta()*rad_to_deg > anal.getMaxTheta(r[1].get3Momentum().Mag(), 1) ) i = 0;
			else if ( rho_mode==0 && !isGoodRho[0] && r[0].get3Momentum().Theta()*rad_to_deg < anal.getMinTheta(r[0].get3Momentum().Mag(), 2) ) i = 1;
			else { continue; }
			

			int chargeIdx = (int)( r[i].getCharge() < 1 );
			double p_pi = r[i].get3Momentum().Mag();

			if( !applyMatching(r[i], p_pi) ) continue;

			if( (*Mx_2pi)*(*Mx_2pi) > Mx2_cut ) continue;
			const int index = 2;

			auto b = calcBins(*e_r, r[i], p_pi);
			if (b.E < 0 || b.Q2 < 0 || b.xB < 0 || b.Z < 0 || b.var < 0) continue;
			corrector.setKinematics( e_r->getXb(), e_r->getQ2(), r[i].getZ(), p_pi );
			if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
			double weight = rhoWeight[i];
			if( (applyCorr & 3) == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( (applyCorr & 3) == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( (applyCorr & 3) == 3 ) weight *= corrector.getCorrectionFactor(2, chargeIdx);
			//if( applyCorr > 3 && b.p >= 0) weight *= corrector.getCorrectionFactor(3, chargeIdx);

			int Q2_lo = (overflowBinning & 1) ? 0     : b.Q2;
			int Q2_hi = b.Q2 + 1;
			int xB_lo = (overflowBinning & 2) ? 0     : b.xB;
			int xB_hi = b.xB + 1;
			int Z_lo  = (overflowBinning & 4) ? 0     : b.Z;
			int Z_hi  = b.Z  + 1;
			for (int iq = Q2_lo; iq < Q2_hi; iq++) {
				for (int ix = xB_lo; ix < xB_hi; ix++) {
					for (int iz = Z_lo; iz < Z_hi; iz++) {
						events_in_bin[index][chargeIdx][b.var][iq][ix][iz][b.p][b.E] += rhoWeight[i];
						weights_in_bin[index][chargeIdx][b.var][iq][ix][iz][b.E] += weight;
						errors_in_bin[index][chargeIdx][b.var][iq][ix][iz][b.p][b.E] += pow(rhoWeight[i],2) + pow(rhoError[i],2);
						
						//events_in_bin[index][1][b.var][iq][ix][iz][b.p][b.E] += rhoWeight[i];
						//weights_in_bin[index][1][b.var][iq][ix][iz][b.E] += weight;
						//errors_in_bin[index][1][b.var][iq][ix][iz][b.p][b.E] += pow(rhoWeight[i],2) + pow(rhoError[i],2);

						// Var(W^h) = sum_i (delta_alpha_i^h)^2  — formula (2)
						// (Pure weight-uncertainty piece, no Poisson on alpha_i.)
						vm_weight_var_in_bin[chargeIdx][b.var][iq][ix][iz][b.E] += pow(rhoError[i], 2);
						//vm_weight_var_in_bin[1][b.var][iq][ix][iz][b.E] += pow(rhoError[i], 2);
					}
				}
			}
		}
	}




	//////////////////////////////////////
	///////////// Calcs //////////////////
	//////////////////////////////////////
	//
	// Full analytical statistical uncertainty on the fragmentation
	// function FF = (4 - R) / (4R - 1), where R = sigma+/sigma-.
	//
	// Per-bin cross-section definition (per charge h in {+,-}):
	//
	//     sigma^h = C_mc^h * sum_p [ C_pi2k^h(p) * N_pi^h(p)
	//                              - C_k2pi^h(p) * N_K^h(p) ]
	//             - sum_{i in n1} C_mc^h * alpha_i^h
	//             - sum_{i in n2} C_mc^{-h} * alpha_i^{-h}
	//
	// where the cross-charge VM term keeps the opposite C_mc factor.
	// Define:
	//     S_PID^h = sum_p [ C_pi2k^h(p)*N_pi^h(p)  -  C_k2pi^h(p)*N_K^h(p) ]
	//     W^h     = sum_i alpha_i^h
	//     B^h     = S_PID^h - W^h
	//
	// Seven-quantity decomposition implemented below:
	//   (1) Var(S_PID^h) — PID + Poisson, per charge
	//   (2) W^h, Var(W^h) — VM weight sums and variances, per charge
	//   (3) B^h, Var(B^h) = Var(S_PID^h) + Var(W^h)
	//   (4) Var(sigma+), Var(sigma-) — cross-section variances
	//   (5) Cov(sigma+, sigma-) — correlation via shared VM weights
	//   (6) Var(R) = R^2 * [ Var(sigma+)/sigma+^2 + Var(sigma-)/sigma-^2
	//                      - 2*Cov/(sigma+ * sigma-) ]
	//   (7) Var(FF) = ( 17 / (4R - 1)^2 )^2 * Var(R)
	//
	// Assumptions:
	//   • PID corrections C_pi2k^h(p) and C_k2pi^h(p) are independent
	//     across p-bins AND across charges (+/-).
	//   • VM weights alpha_i^h are independent across events AND across
	//     charges; within a charge, their sum W^h is treated as a sum
	//     of independent random variables.
	//   • C_mc+ and C_mc- are independent (different MC samples or
	//     independent binwise fluctuations).
	//   • Raw counts N_pi^h(p), N_K^h(p) follow Poisson statistics:
	//     Var(N) = <N>, so the counting contribution to Var(S_PID^h)
	//     is (C_pi2k)^2 * N_pi + (C_k2pi)^2 * N_K.
	//   • The only source of correlation between sigma+ and sigma- is
	//     the cross-charge VM term sharing W^{-h} and the same C_mc+/-.
	//

	cout << "Start calculations\n";

	// ── RICH diagnostic summary ──
	{
		long long total_pi = 0, rich_pi = 0, total_K = 0, rich_K = 0;
		for (int c = 0; c < 2; c++)
			for (int v = 0; v < bins_var; v++)
				for (int q = 0; q < bins_Q2; q++)
					for (int x = 0; x < bins_xB; x++)
						for (int z = 0; z < bins_Z; z++)
							for (int p = 0; p < bins_p; p++)
								for (int e = 0; e < N_E; e++) {
									total_pi += (long long)events_in_bin[0][c][v][q][x][z][p][e];
									rich_pi  += (long long)rich_in_bin[0][c][v][q][x][z][p][e];
									total_K  += (long long)events_in_bin[1][c][v][q][x][z][p][e];
									rich_K   += (long long)rich_in_bin[1][c][v][q][x][z][p][e];
								}
		cout << "RICH diagnostic: pions with RICH = " << rich_pi << " / " << total_pi
			 << " (" << (total_pi > 0 ? 100.*rich_pi/total_pi : 0) << "%)\n";
		cout << "RICH diagnostic: kaons with RICH = " << rich_K << " / " << total_K
			 << " (" << (total_K > 0 ? 100.*rich_K/total_K : 0) << "%)\n";
		if (rich_pi == 0 && total_pi > 0)
			cout << "WARNING: No RICH info found for pions. delta-T reverts to pure Poisson (k treated as fixed).\n";
	}

	// ── Diagnostic histograms for uncertainty decomposition ──
	// hSigma[var][Q2][xB][charge]: sigma+/- with delta-sigma as error bar
	// hCovHist[var][Q2][xB]:       Cov(sigma+, sigma-) per z-bin
	// hDeltaT[var][Q2][xB][charge]: 2D (z x p) showing sqrt(delta-T^2) per bin
	TH1F * hSigma[bins_var][bins_Q2][bins_xB][2];
	TH1F * hCovHist[bins_var][bins_Q2][bins_xB];
	TH2F * hDeltaT[bins_var][bins_Q2][bins_xB][2];
	// Fragmentation function FF = (4 - R) / (4R - 1)  (formula (7) context)
	TH1F * hFF[bins_var][bins_Q2][bins_xB];

	for (int var = 0; var < bins_var; var++) {
		for (int ii = 0; ii < bins_Q2; ii++) {
			for (int jj = 0; jj < bins_xB; jj++) {
				TString ds = Form("_%d_%d_%d", var, ii+1, jj+1);
				for (int c = 0; c < 2; c++) {
					hSigma[var][ii][jj][c] = new TH1F(
						Form("hSigma%s%s%s", (c==0 ? "Plus" : "Minus"),
							 charge_str[c].Data(), ds.Data()),
						"", bins_Z, 0.3, 1.0);
					hSigma[var][ii][jj][c]->Sumw2();

					hDeltaT[var][ii][jj][c] = new TH2F(
						Form("hDeltaT%s%s", charge_str[c].Data(), ds.Data()),
						";z;p bin", bins_Z, 0.3, 1.0, bins_p, 0, bins_p);
					hDeltaT[var][ii][jj][c]->Sumw2();
				}
				hCovHist[var][ii][jj] = new TH1F(
					Form("hCov%s", ds.Data()), "", bins_Z, 0.3, 1.0);

				// Fragmentation function per z-bin (one per Q2/xB/var)
				hFF[var][ii][jj] = new TH1F(
					Form("hFF%s", ds.Data()),
					";z;FF = (4 - R)/(4R - 1)",
					bins_Z, 0.3, 1.0);
				hFF[var][ii][jj]->Sumw2();
			}
		}
	}

	// ── Determine which MC correction type to use ──
	// bits[1:0] of applyCorr: 0=none, 1=bin migration, 2=acceptance, 3=full MC
	int mc_factor_type = (applyCorr & 3) - 1;  // -1 = no correction (C_MC = 1)

	double xB_step = (xB_max - xB_min) / (double)bins_xB;
	double Q2_step = (Q2_max - Q2_min) / (double)bins_Q2;

	// ── Determine which beam energies have data ──
	// Only load correction files for energies that were actually filled.
	// (Mirrors the original setBin guard: "if (weights == 0) return;"
	//  which skipped updateCorrectionsForBeam for empty bins.)
	bool beamHasData[N_E] = {false, false, false};
	for (int E = 0; E < N_E; E++)
		for (int c = 0; c < 2 && !beamHasData[E]; c++)
			for (int v = 0; v < bins_var && !beamHasData[E]; v++)
				for (int q = 0; q < bins_Q2 && !beamHasData[E]; q++)
					for (int x = 0; x < bins_xB && !beamHasData[E]; x++)
						for (int z = 0; z < bins_Z && !beamHasData[E]; z++)
							for (int p = 0; p < bins_p && !beamHasData[E]; p++)
								if (events_in_bin[0][c][v][q][x][z][p][E] > 0)
									beamHasData[E] = true;

	for (int E = 0; E < N_E; E++)
		cout << "Beam energy E=" << E << " (Ebeam=" << 10.2 + 0.2*E
			 << "): " << (beamHasData[E] ? "HAS DATA" : "no data, skipping") << "\n";

	// ── Pre-load corrections once per beam energy ──────────────────────────────
	// Loading involves cloning TF3 objects (expensive formula re-parse). Doing it
	// inside the (Q2, xB, var) loop would repeat it O(bins_Q2*bins_xB*bins_var)
	// times; pre-loading here reduces it to at most N_E loads total.
	correctionTools* corrByE[N_E] = {};
	for (int E = 0; E < N_E; E++) {
		if (!beamHasData[E]) continue;
		corrByE[E] = new correctionTools(n4d_corr_bins > 0 ? 4 : 2);
		corrByE[E]->setN4dBins(n4d_corr_bins);
		corrByE[E]->setMisidFitName("corrections_misid_fit.root");
		TString corrName_e = correction_name;
		updateCorrectionsForBeam(E, corrName_e, *corrByE[E]);
		cout << "Pre-loaded corrections for E=" << E
		     << " (Ebeam=" << 10.2 + 0.2*E << ")\n";
	}

	// ── Main calculation loop ──
	for (int i = 1; i <= bins_Q2; i++) {
		for (int j = 1; j <= bins_xB; j++) {
			for (int var = 0; var < bins_var; var++) {

				double xB_center = xB_min + (j - 0.5) * xB_step;
				double Q2_center = Q2_min + (i - 0.5) * Q2_step;

				// Per-z-bin accumulators (summed over beam energies E)
				double sigma_total[2][bins_Z];        // sigma+ and sigma-
				double dsigma_sq_total[2][bins_Z];    // (delta-sigma+-)^2
				double cov_total_z[bins_Z];           // Cov(sigma+, sigma-)
				double dT_sq_accum[2][bins_Z][bins_p]; // per-bin delta-T^2 diagnostic

				// Per-term breakdown of Var(sigma^h) for each charge [c][t][z].
				// See Formula (4) block below for the full definition; indices:
				//   [0] = (B^h)^2 * sigma_{C_mc^h}^2        (C_mc, same-charge)
				//   [1] = (C_mc^h)^2 * Var(S_PID^h)         (counting + PID-corr)
				//   [2] = (C_mc^h)^2 * Var(W^h)             (VM, same-charge)
				//   [3] = (W^{-h})^2 * sigma_{C_mc^{-h}}^2  (C_mc, cross-charge)
				//   [4] = (C_mc^{-h})^2 * Var(W^{-h})       (VM, cross-charge)
				double dsig_term[2][5][bins_Z];
				// Categorized breakdown of Var(sigma^h) [c][category][z]
				// exposing the four sources the user requested:
				//   [0] counting-statistics : (C_mc^h)^2 * Var_S_PID_poisson^h
				//   [1] PID-correction      : (C_mc^h)^2 * Var_S_PID_pidcorr^h
				//   [2] VM-weight           : (C_mc^h)^2 * Var(W^h)
				//                           + (C_mc^{-h})^2 * Var(W^{-h})
				//   [3] C_mc                : (B^h)^2 * sigma_{C_mc^h}^2
				//                           + (W^{-h})^2 * sigma_{C_mc^{-h}}^2
				double dsig_cat[2][4][bins_Z];
				// Covariance term breakdown [t][z] (see Formula (5)):
				//   [0] = -B+ * W+ * sigma_{C_mc+}^2
				//   [1] = -B- * W- * sigma_{C_mc-}^2
				//   [2] = (C_mc+)^2 * Var(W+)
				//   [3] = (C_mc-)^2 * Var(W-)
				double cov_term[4][bins_Z];

				memset(sigma_total,     0, sizeof(sigma_total));
				memset(dsigma_sq_total, 0, sizeof(dsigma_sq_total));
				memset(cov_total_z,     0, sizeof(cov_total_z));
				memset(dT_sq_accum,     0, sizeof(dT_sq_accum));
				memset(dsig_term,       0, sizeof(dsig_term));
				memset(dsig_cat,        0, sizeof(dsig_cat));
				memset(cov_term,        0, sizeof(cov_term));

				for (int E = 0; E < N_E; E++) {
					if (!beamHasData[E] || !corrByE[E]) continue;
					correctionTools& corrE = *corrByE[E];

					for (int zbin = 0; zbin < bins_Z; zbin++) {
						double z_center = hZ[var][i-1][j-1][0]->GetXaxis()->GetBinCenter(zbin+1);

						// Use event-weighted mean kinematics when overflow
						// binning is active; fall back to bin center otherwise.
						double cnt = count_kin[var][i-1][j-1][zbin][E];
						double Q2_kin = (cnt > 0) ? sum_Q2_kin[var][i-1][j-1][zbin][E] / cnt : Q2_center;
						double xB_kin = (cnt > 0) ? sum_xB_kin[var][i-1][j-1][zbin][E] / cnt : xB_center;
						double z_kin  = (cnt > 0) ? sum_z_kin [var][i-1][j-1][zbin][E] / cnt : z_center;

						// ── Get C_MC+- and delta-C_MC+- ──
						corrE.setKinematics(xB_kin, Q2_kin, z_kin, 1.5);
						if (n4d_corr_bins > 0) corrE.set4dBin(var);

						double C_MC[2], dC_MC[2];
						for (int c = 0; c < 2; c++) {
							C_MC[c]  = (mc_factor_type >= 0) ?
								corrE.getCorrectionFactor(mc_factor_type, c) : 1.0;
							dC_MC[c] = (mc_factor_type >= 0) ?
								corrE.getCorrectionError(mc_factor_type, c) : 0.0;
						}

						// =====================================================
						// Formula (1):  S_PID^h and Var(S_PID^h), per charge
						//
						//   S_PID^h = sum_p [ C_pi2k^h(p)*N_pi^h(p)
						//                   - C_k2pi^h(p)*N_K^h(p) ]
						//
						//   Var(S_PID^h) = sum_p [
						//         (N_pi^h(p))^2 * sigma_{C_pi2k^h(p)}^2
						//       + (C_pi2k^h(p))^2 * N_pi^h(p)
						//       + (N_K^h(p))^2  * sigma_{C_k2pi^h(p)}^2
						//       + (C_k2pi^h(p))^2 * N_K^h(p) ]
						//
						// (PID corrections independent across p and across
						//  charges; Poisson counting on N_pi, N_K.)
						// =====================================================
						double S_PID[2]             = {0, 0};
						double Var_S_PID[2]         = {0, 0};
						// Split of Var(S_PID^h) for diagnostic output:
						//   _pidcorr : (N_pi)^2 * sigma_{C_pi2k}^2 + (N_K)^2 * sigma_{C_k2pi}^2
						//   _poisson : (C_pi2k)^2 * N_pi           + (C_k2pi)^2 * N_K
						double Var_S_PID_pidcorr[2] = {0, 0};
						double Var_S_PID_poisson[2] = {0, 0};

						for (int pbin = 0; pbin < bins_p; pbin++) {
							double p_center = (p_bin_edges[pbin] + p_bin_edges[pbin+1]) / 2.;
							corrE.setKinematics(xB_kin, Q2_kin, z_kin, p_center);
							if (n4d_corr_bins > 0) corrE.set4dBin(var);

							for (int c = 0; c < 2; c++) {
								double N_pi = events_in_bin[0][c][var][i-1][j-1][zbin][pbin][E];
								double N_K  = (applyCorr & 4) ?
									events_in_bin[1][c][var][i-1][j-1][zbin][pbin][E] : 0.0;

								// C_pi2k^h(p) and its uncertainty (pi RICH purity)
								double C_pi2k  = (applyCorr & 4) ?
									corrE.getCorrectionFactor(3, c) : 1.0;
								double dC_pi2k = (applyCorr & 4) ?
									corrE.getCorrectionError(3, c) : 0.0;
								// C_k2pi^h(p) and its uncertainty (K->pi misid subtraction)
								double C_k2pi  = (applyCorr & 4) ?
									corrE.getCorrectionFactor(4, c) : 0.0;
								double dC_k2pi = (applyCorr & 4) ?
									corrE.getCorrectionError(4, c) : 0.0;

								// S_PID^h contribution at this p-bin
								// (sign: +C_pi2k*N_pi recovers true pions,
								//        -C_k2pi*N_K  subtracts K->pi misid)
								S_PID[c] += C_pi2k * N_pi - C_k2pi * N_K;

								// Per-p variance contribution to Var(S_PID^h),
								// split into PID-correction and Poisson pieces:
								double v_pidcorr = N_pi * N_pi * dC_pi2k * dC_pi2k
								                 + N_K  * N_K  * dC_k2pi * dC_k2pi;
								double v_poisson = C_pi2k * C_pi2k * N_pi
								                 + C_k2pi * C_k2pi * N_K;
								double v_p = v_pidcorr + v_poisson;

								Var_S_PID[c]               += v_p;
								Var_S_PID_pidcorr[c]       += v_pidcorr;
								Var_S_PID_poisson[c]       += v_poisson;
								dT_sq_accum[c][zbin][pbin] += v_p;
							} // charge
						} // pbin

						// =====================================================
						// Formula (2):  W^h and Var(W^h), per charge
						//
						//   W^h      = sum_i alpha_i^h
						//   Var(W^h) = sum_i (delta_alpha_i^h)^2
						//
						// W^h is already accumulated at fill time in
						// events_in_bin[2] (summed over p); Var(W^h) is
						// accumulated separately in vm_weight_var_in_bin
						// (pure weight-uncertainty piece, no Poisson on alpha_i).
						// =====================================================
						double W[2]     = {0, 0};
						double Var_W[2] = {0, 0};

						if (applyCorr & 8) {
							for (int c = 0; c < 2; c++) {
								for (int pbin = 0; pbin < bins_p; pbin++) {
									W[c] += events_in_bin[2][c][var][i-1][j-1][zbin][pbin][E];
								}
								Var_W[c] += vm_weight_var_in_bin[c][var][i-1][j-1][zbin][E];
							}
						}

						// =====================================================
						// Formula (3):  B^h = S_PID^h - W^h, Var(B^h)
						//
						//   Var(B^h) = Var(S_PID^h) + Var(W^h)
						//
						// (S_PID^h uses per-event RICH-tagged counts, W^h uses
						//  per-event VM weights — assumed independent.)
						// =====================================================
						double B[2], Var_B[2];
						for (int c = 0; c < 2; c++) {
							B[c]     = S_PID[c] - W[c];
							Var_B[c] = Var_S_PID[c] + Var_W[c];
						}
						// Var_B enters Formula (4) as (C_mc^h)^2 * Var(B^h), which
						// we decompose below into its Var(S_PID^h) + Var(W^h) pieces
						// for diagnostic granularity.  Reference here to document
						// the identity and silence unused-variable warnings.
						(void)Var_B[0]; (void)Var_B[1];

						// =====================================================
						// Formula (4):  Var(sigma+) and Var(sigma-)
						//
						//   sigma^h = C_mc^h * B^h  -  C_mc^{-h} * W^{-h}
						//
						//   Var(sigma+) = (C_mc+)^2 * Var(B+) + (B+)^2 * sigma_{C_mc+}^2
						//               + (C_mc-)^2 * Var(W-) + (W-)^2 * sigma_{C_mc-}^2
						//   Var(sigma-) : swap +<->-
						//
						// Per-term breakdown indices (stored in dsig_term[c][t][z]):
						//   [0] = (B^h)^2 * sigma_{C_mc^h}^2        (C_mc, same-charge)
						//   [1] = (C_mc^h)^2 * Var(S_PID^h)         (counting stats)
						//   [2] = (C_mc^h)^2 * Var(W^h)             (VM, same-charge)
						//   [3] = (W^{-h})^2 * sigma_{C_mc^{-h}}^2  (C_mc, cross-charge)
						//   [4] = (C_mc^{-h})^2 * Var(W^{-h})       (VM, cross-charge)
						// Sum of these equals Var(sigma^h), with Var(B^h) split
						// into its Var(S_PID^h) and Var(W^h) pieces for diagnosis.
						// =====================================================
						double sig[2];
						sig[0] = C_MC[0] * B[0] - C_MC[1] * W[1];
						sig[1] = C_MC[1] * B[1] - C_MC[0] * W[0];

						double t0[2][5];
						// sigma+ terms (c = 0)
						t0[0][0] = pow(B[0], 2)    * pow(dC_MC[0], 2);
						t0[0][1] = pow(C_MC[0], 2) * Var_S_PID[0];
						t0[0][2] = pow(C_MC[0], 2) * Var_W[0];
						t0[0][3] = pow(W[1], 2)    * pow(dC_MC[1], 2);
						t0[0][4] = pow(C_MC[1], 2) * Var_W[1];
						// sigma- terms (c = 1)
						t0[1][0] = pow(B[1], 2)    * pow(dC_MC[1], 2);
						t0[1][1] = pow(C_MC[1], 2) * Var_S_PID[1];
						t0[1][2] = pow(C_MC[1], 2) * Var_W[1];
						t0[1][3] = pow(W[0], 2)    * pow(dC_MC[0], 2);
						t0[1][4] = pow(C_MC[0], 2) * Var_W[0];

						double dsig_sq[2] = {0, 0};
						for (int c = 0; c < 2; c++)
							for (int t = 0; t < 5; t++) {
								dsig_sq[c]            += t0[c][t];
								dsig_term[c][t][zbin] += t0[c][t];
							}

						// Four-source categorized breakdown (see dsig_cat header):
						//   [0] counting stats  : (C_mc^h)^2 * Var_S_PID_poisson^h
						//   [1] PID correction  : (C_mc^h)^2 * Var_S_PID_pidcorr^h
						//   [2] VM weight       : same-charge + cross-charge VM terms
						//   [3] C_mc            : same-charge + cross-charge C_mc terms
						for (int c = 0; c < 2; c++) {
							int oc = 1 - c;  // opposite charge index
							dsig_cat[c][0][zbin] += pow(C_MC[c], 2) * Var_S_PID_poisson[c];
							dsig_cat[c][1][zbin] += pow(C_MC[c], 2) * Var_S_PID_pidcorr[c];
							dsig_cat[c][2][zbin] += pow(C_MC[c],  2) * Var_W[c]
							                      + pow(C_MC[oc], 2) * Var_W[oc];
							dsig_cat[c][3][zbin] += pow(B[c],  2) * pow(dC_MC[c],  2)
							                      + pow(W[oc], 2) * pow(dC_MC[oc], 2);
						}

						// =====================================================
						// Formula (5):  Cov(sigma+, sigma-)
						//
						// Correlation sources: shared C_mc+ (couples the B+
						// term in sigma+ to the W+ term in sigma- via dC_mc+),
						// shared C_mc-, and shared W^h appearing in both
						// sigma^h (inside B^h) and sigma^{-h} (as cross-charge
						// subtraction).
						//
						//   Cov(sigma+, sigma-) = -B+ * W+ * sigma_{C_mc+}^2
						//                        -B- * W- * sigma_{C_mc-}^2
						//                        +(C_mc+)^2 * Var(W+)
						//                        +(C_mc-)^2 * Var(W-)
						// =====================================================
						double ct[4];
						ct[0] = -B[0] * W[0] * pow(dC_MC[0], 2);
						ct[1] = -B[1] * W[1] * pow(dC_MC[1], 2);
						ct[2] =  pow(C_MC[0], 2) * Var_W[0];
						ct[3] =  pow(C_MC[1], 2) * Var_W[1];

						double cov_E = ct[0] + ct[1] + ct[2] + ct[3];

						// ── Accumulate across beam energies ──
						// (independent beam energies => variances add)
						for (int c = 0; c < 2; c++) {
							sigma_total[c][zbin]      += sig[c];
							dsigma_sq_total[c][zbin]  += dsig_sq[c];
						}
						cov_total_z[zbin] += cov_E;
						for (int t = 0; t < 4; t++)
							cov_term[t][zbin] += ct[t];

					} // zbin
				} // E

				// ══════════════════════════════════════════════════════════════
				// SET BIN CONTENT — output histograms filled here
				//
				// Inputs come from the E-loop above (sigma_total, dsigma_sq_total,
				// cov_total_z), which accumulated contributions from each beam energy.
				// The following block translates those per-z accumulators into the
				// final histogram bin values written to the output file.
				// ══════════════════════════════════════════════════════════════
				for (int zbin = 0; zbin < bins_Z; zbin++) {
					double sp     = sigma_total[0][zbin];
					double sm     = sigma_total[1][zbin];
					double dsp_sq = dsigma_sq_total[0][zbin];
					double dsm_sq = dsigma_sq_total[1][zbin];
					double cov    = cov_total_z[zbin];

					// ── Bin content: sigma+/- (cross-sections, with statistical error) ──
					hSigma[var][i-1][j-1][0]->SetBinContent(zbin+1, sp);
					hSigma[var][i-1][j-1][0]->SetBinError(zbin+1,
						sqrt(fmax(dsp_sq, 0.0)));
					hSigma[var][i-1][j-1][1]->SetBinContent(zbin+1, sm);
					hSigma[var][i-1][j-1][1]->SetBinError(zbin+1,
						sqrt(fmax(dsm_sq, 0.0)));
					// ── Bin content: Cov(sigma+, sigma-) ──
					hCovHist[var][i-1][j-1]->SetBinContent(zbin+1, cov);

					// =====================================================
					// Formula (6):  Var(R), R = sigma+ / sigma-
					//
					//   Var(R) = R^2 * [ Var(sigma+)/sigma+^2
					//                  + Var(sigma-)/sigma-^2
					//                  - 2*Cov(sigma+, sigma-) / (sigma+*sigma-) ]
					//
					// Formula (7):  Var(FF), FF = (4 - R) / (4R - 1)
					//
					//   dFF/dR = -17 / (4R - 1)^2
					//   Var(FF) = ( 17 / (4R - 1)^2 )^2 * Var(R)
					//
					// Guards: sigma+, sigma-, (4R - 1) must be nonzero and
					// finite; otherwise the bin is flagged as empty/pathological
					// and bin content/error set to 0.
					// =====================================================
					const double kEps = 1e-12;
					bool sigma_ok = (isfinite(sp) && isfinite(sm)
					              && fabs(sp) > kEps && fabs(sm) > kEps);

					if (sigma_ok) {
						double R = sp / sm;

						// ── Formula (6): Var(R) ──
						double dR_over_R_sq = dsp_sq / (sp * sp)
						                    + dsm_sq / (sm * sm)
						                    - 2.0 * cov / (sp * sm);
						double sigma_R_sq = R * R * fmax(dR_over_R_sq, 0.0);
						double dR = sqrt(sigma_R_sq);

						// ── Bin content: R = sigma+/sigma- (multiplicity ratio) ──
						hZ[var][i-1][j-1][0]->SetBinContent(zbin+1, R);
						hZ[var][i-1][j-1][0]->SetBinError(zbin+1, dR);

						// ── Formula (7): Var(FF) ──
						double denom = 4.0 * R - 1.0;
						if (isfinite(R) && fabs(denom) > kEps) {
							double FF = (4.0 - R) / denom;
							if (isfinite(FF)) {
								double dFF = (17.0 / (denom * denom)) * sqrt(fmax(sigma_R_sq, 0.0));
								// ── Bin content: FF = (4 - R)/(4R - 1) (fragmentation function) ──
								hFF[var][i-1][j-1]->SetBinContent(zbin+1, FF);
								hFF[var][i-1][j-1]->SetBinError  (zbin+1, dFF);
							} else {
								hFF[var][i-1][j-1]->SetBinContent(zbin+1, 0);
								hFF[var][i-1][j-1]->SetBinError  (zbin+1, 0);
							}
						} else {
							// (4R - 1) ~ 0 or R not finite: FF undefined
							hFF[var][i-1][j-1]->SetBinContent(zbin+1, 0);
							hFF[var][i-1][j-1]->SetBinError  (zbin+1, 0);
						}
					} else {
						// sigma+ or sigma- is zero/non-finite: set R and FF to 0
						hZ [var][i-1][j-1][0]->SetBinContent(zbin+1, 0);
						hZ [var][i-1][j-1][0]->SetBinError  (zbin+1, 0);
						hFF[var][i-1][j-1]   ->SetBinContent(zbin+1, 0);
						hFF[var][i-1][j-1]   ->SetBinError  (zbin+1, 0);
					}

					// ── Bin content: delta-T diagnostic (sqrt of summed variance over E, per z×p cell) ──
					for (int c = 0; c < 2; c++)
						for (int pbin = 0; pbin < bins_p; pbin++)
							hDeltaT[var][i-1][j-1][c]->SetBinContent(
								hDeltaT[var][i-1][j-1][c]->GetBin(zbin+1, pbin+1),
								sqrt(fmax(dT_sq_accum[c][zbin][pbin], 0.0)));
				} // zbin

				// ── Diagnostic: fractional contribution of each variance source ──
				// Gated behind the `verbose` flag (CLI arg 18).  Prints, per bin,
				// the fractional contribution to Var(sigma+) and Var(sigma-) of
				// the four sources: counting statistics, PID correction,
				// VM weights, and C_mc; plus the covariance and (dR/R)^2
				// breakdowns for completeness.
				if (verbose) {
					// Accumulate term sums across z-bins
					double sum_dsig[2]     = {0, 0};
					double sum_term[2][5]  = {{0}};
					double sum_cat[2][4]   = {{0}};
					double sum_ct[4]       = {0};
					double sum_dRoR_sq     = 0;
					double sum_dRoR_sp = 0, sum_dRoR_sm = 0, sum_dRoR_cov = 0;
					int nz_filled = 0;

					for (int zbin = 0; zbin < bins_Z; zbin++) {
						double sp = sigma_total[0][zbin];
						double sm = sigma_total[1][zbin];
						if (sp == 0 || sm == 0) continue;
						nz_filled++;

						for (int c = 0; c < 2; c++) {
							sum_dsig[c] += dsigma_sq_total[c][zbin];
							for (int t = 0; t < 5; t++)
								sum_term[c][t] += dsig_term[c][t][zbin];
							for (int k = 0; k < 4; k++)
								sum_cat[c][k] += dsig_cat[c][k][zbin];
						}
						double cov = cov_total_z[zbin];
						for (int t = 0; t < 4; t++)
							sum_ct[t] += cov_term[t][zbin];

						// Three master-formula contributions
						double dsp_sq = dsigma_sq_total[0][zbin];
						double dsm_sq = dsigma_sq_total[1][zbin];
						sum_dRoR_sp  += dsp_sq / (sp * sp);
						sum_dRoR_sm  += dsm_sq / (sm * sm);
						sum_dRoR_cov += -2.0 * cov / (sp * sm);
						sum_dRoR_sq  += dsp_sq / (sp * sp)
						              + dsm_sq / (sm * sm)
						              - 2.0 * cov / (sp * sm);
					}

					if (nz_filled > 0) {
						const char* catname[4] = {
							"counting stats  ",
							"PID correction  ",
							"VM weights      ",
							"C_mc correction "
						};
						const char* tname[5] = {
							"C_mc same-chg   ", "counting+PID-corr",
							"VM same-chg     ", "C_mc cross-chg  ",
							"VM cross-chg    "
						};
						const char* cname[4] = {
							"Cov: -B+  * W+  * dC+^2",
							"Cov: -B-  * W-  * dC-^2",
							"Cov: +C+^2 * Var(W+)   ",
							"Cov: +C-^2 * Var(W-)   "
						};
						cout << "\n──── Uncertainty breakdown: Q2 bin " << i
						     << ", xB bin " << j << ", var " << var
						     << " (" << nz_filled << " z-bins) ────\n";

						// Four-source breakdown (user-requested)
						for (int c = 0; c < 2; c++) {
							cout << "  Var(sigma" << (c==0 ? "+" : "-")
							     << ") — 4-source breakdown:\n";
							if (sum_dsig[c] > 0)
								for (int k = 0; k < 4; k++)
									cout << "    " << catname[k] << ":  "
									     << Form("%6.2f%%", 100. * sum_cat[c][k] / sum_dsig[c])
									     << "\n";
						}

						// Full 5-term breakdown (for deeper debugging)
						for (int c = 0; c < 2; c++) {
							cout << "  Var(sigma" << (c==0 ? "+" : "-")
							     << ") — 5-term detail:\n";
							if (sum_dsig[c] > 0)
								for (int t = 0; t < 5; t++)
									cout << "    " << tname[t] << ":  "
									     << Form("%6.2f%%", 100. * sum_term[c][t] / sum_dsig[c])
									     << "\n";
						}

						cout << "  Cov(sigma+, sigma-) breakdown:\n";
						double sum_ct_abs = fabs(sum_ct[0]) + fabs(sum_ct[1])
						                  + fabs(sum_ct[2]) + fabs(sum_ct[3]);
						if (sum_ct_abs > 0)
							for (int t = 0; t < 4; t++)
								cout << "    " << cname[t] << ":  "
								     << Form("%+7.2f%%", 100. * sum_ct[t] / sum_ct_abs)
								     << " (of |Cov| total)\n";

						cout << "  (dR/R)^2 breakdown (summed over z):\n";
						if (sum_dRoR_sq > 0) {
							cout << "    (dsig+/sig+)^2 :  "
							     << Form("%6.2f%%", 100. * sum_dRoR_sp / sum_dRoR_sq) << "\n";
							cout << "    (dsig-/sig-)^2 :  "
							     << Form("%6.2f%%", 100. * sum_dRoR_sm / sum_dRoR_sq) << "\n";
							cout << "    -2*Cov/(s+*s-) :  "
							     << Form("%+6.2f%%", 100. * sum_dRoR_cov / sum_dRoR_sq) << "\n";
						}
						cout << "\n";
					}
				}

				// Write all histograms for this kinematic bin
				outFile->cd();
				hZ[var][i-1][j-1][0]->Write();
				hSigma[var][i-1][j-1][0]->Write();
				hSigma[var][i-1][j-1][1]->Write();
				hCovHist[var][i-1][j-1]->Write();
				hDeltaT[var][i-1][j-1][0]->Write();
				hDeltaT[var][i-1][j-1][1]->Write();
				hFF[var][i-1][j-1]->Write();   // Fragmentation function (Formula (7))

			} // var
		} // j (xB)
	} // i (Q2)
	/*
	TH1F * rat = new TH1F("dOu", "dOu", bins_xB, xB_min, xB_max);
	for( int x = 0; x < bins_xB; x++ ){
		for( int q = 4; q<bins_Q2; q++ ){
			hZ[0][3][x][0]->Add(hZ[0][q][x][0]);
			hZ[0][3][x][1]->Add(hZ[0][q][x][1]);
		}
		double rat_bin_p = 0;
		double rat_bin_m = 0;
		for( int k = 2; k <= bins_Z; k++ ){
			if( isnan(hZ[0][3][x][0]->GetBinContent(k)) || isnan(hZ[0][3][x][1]->GetBinContent(k))) continue;
			rat_bin_p+=hZ[0][3][x][0]->GetBinContent(k);
			rat_bin_m+=hZ[0][3][x][1]->GetBinContent(k);
		}
		cout<<"pip: "<<rat_bin_p<<std::endl;
		cout<<"pim: "<<rat_bin_m<<std::endl;

		rat->SetBinContent(x + 1, rat_bin_p/rat_bin_m);
		
		
	}
	*/
	//rat->Write();

	
	for (int E = 0; E < N_E; E++) delete corrByE[E];

	cout<<"Closing Out File\n";
	delete outFile;
	cout << "Done\n";

}


void zeroSuppress( TH1F * h){

	int nBins = h->GetNbinsX();
	for( int i = 1; i <= nBins; i++ ){
		if( h->GetBinContent(i) < 0 ){
			h->SetBinContent(i, 0);
			h->SetBinError(i, 0);
		}
	}
}


double getVarVal(TString var,  electron  e, pion pi ){
	if( var == "null" ) return 0;
	if( var == "p_e" ) return e.get3Momentum().Mag();
	if( var == "theta_e" ) return e.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e.get3Momentum().Phi()*rad_to_deg;
	if( var == "W2" ) return e.getW2();
	if( var == "Q2" ) return e.getQ2();
	if( var == "xB" ) return e.getXb();
	if( var == "y" ) return e.getY();
	if( var == "sector" ) return e.getDC_sector();
	if( var == "sector_e" ) { double phi_deg = e.get3Momentum().Phi()*rad_to_deg + 20.; if (phi_deg < 0.) phi_deg += 360.; return (int)(phi_deg / 60.) + 1; }


	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) {
		return pi.getPi_q().Phi()*rad_to_deg + 360*( (int) (pi.getPi_q().Phi() < 0 ) ) ;
	}
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return  pi.getPi_q().Pt();
	if( var == "sector_pi" ) return pi.getDC_sector();
	if( var == "eta" ) return pi.getEta();

	return 0;
}

void setBin(
    TH1F* h,
    int z_bin,
    const VecZ& events,      // [Z][p][E]
    const VecZ& errors,      // [Z][p][E]
    double weights,
    correctionTools& corr,
	int corrType,
    int chargeIdx,
    int E_bin,
    TString corrName
){

	if ( weights == 0 )return;

	double error = 0;
	double term_1 = 0;
	double term_2 = 0;
	double term_3 = 0;
	int E = E_bin;

	double mc_weight = 1;//corr.getCorrectionFactor( 2, chargeIdx );
	double k_weight = 1;

	double mc_err = 0;//corr.getCorrectionError(0, chargeIdx);
	double k_err = 0;//corr.getCorrectionError(1, chargeIdx); 


	updateCorrectionsForBeam(E, corrName, corr);

	//corr.testHists();

	switch ( corrType ){
		case 1:
			mc_weight = corr.getCorrectionFactor(0, chargeIdx);
			mc_err = corr.getCorrectionError(0, chargeIdx);
			break;
		case 2:
			mc_weight = corr.getCorrectionFactor(1, chargeIdx);
			mc_err = corr.getCorrectionError(1, chargeIdx);
			break;
		case 3: case 4: case 5:
			mc_weight = corr.getCorrectionFactor(2, chargeIdx);
			mc_err = corr.getCorrectionError(2, chargeIdx);
			break;
	}
	

	double term_0 = 0; //temp term
	for( int i = 0; i < 4; i++ ){
		double p = (p_bin_edges[i] + p_bin_edges[i+1])/2.;
		corr.setKinematics( corr.getX(), 
										corr.getQ2(),
										corr.getZ(), p);
		switch ( corrType ){
			case 4:
				k_weight = corr.getCorrectionFactor(3, chargeIdx);
				k_err = corr.getCorrectionError(3, chargeIdx);
				break;
			case 5:
				k_weight = corr.getCorrectionFactor(4, chargeIdx);
				k_err = corr.getCorrectionError(4, chargeIdx);
				break;
		
		}
		term_1 += pow( mc_weight*k_weight, 2)*errors[z_bin][i][E];
		term_2 += pow( mc_weight*events[z_bin][i][E]*k_err, 2 );
		term_0 += k_weight*events[z_bin][i][E];
		
		
	}
	term_3 += pow(term_0*mc_err, 2);
	

	error = sqrt( term_1 + term_2 + term_3);

	double bin_old = h->GetBinContent(z_bin+1);
	double err_old = h->GetBinError(z_bin+1);

	h->SetBinContent( z_bin +1, bin_old + weights);
	h->SetBinError(z_bin + 1, sqrt( pow(err_old, 2)+ pow(error,2) ));
}



bool updateCorrectionsForBeam(
    double eBeam,
    double& beam_energy,
    TString& corrName,
    correctionTools& corrector
){
    if (eBeam == beam_energy)
        return false;

	getNewRhoName(eBeam, corrName);
	corrector.setWeightName( corrName);
	if( n4d_corr_bins > 0 ){
		corrector.setWeightName4D( corrName);
		corrector.setN4dBins(n4d_corr_bins);
	}

	// Derive the continuous fit filename from the base template, replacing energy.
	TString fitName = corrName;
	getNewRhoName(eBeam, fitName);
	// Only reload if the fit file name actually changed — avoids re-parsing
	// TF3 formulas (expensive) when the corrections don't depend on beam energy.
	if( fitName != corrector.getWeightFitName() ){
		corrector.setWeightFitName( fitName );
		corrector.loadContinuousCorrections();
	}

    beam_energy = eBeam;
    return true;
}

bool updateCorrectionsForBeam(
	int E,
    TString& corrName,
    correctionTools& corrector
){
	double eBeam = 10.2 + 0.2*E;
	getNewRhoName(eBeam, corrName);
	corrector.setWeightName( corrName);
	if( n4d_corr_bins > 0 ){
		corrector.setWeightName4D( corrName);
		corrector.setN4dBins(n4d_corr_bins);
	}

	// Derive the continuous fit filename from the base template, replacing energy.
	TString fitName = corrName;//"corrections_10.2_fit.root";
	getNewRhoName(eBeam, fitName);
	corrector.setWeightFitName( fitName );
	corrector.loadContinuousCorrections();

    return true;
}

// Helper function to load RhoNorms for multiple beam energies


void getNewRhoName( double E, TString& rho_norm_name) {
    std::string fname = rho_norm_name.Data();  // Convert to std::string for regex
    std::ostringstream ss;
    ss << E;                       // Convert double to string (e.g., 10.6)
    std::string newEnergyStr = ss.str();

    // Regex to match the first occurrence of 10.x or 10.xx
    std::regex energyPattern(R"(10(\.\d+)?)"); 

    // Replace the first match with newEnergyStr
    fname = std::regex_replace(fname, energyPattern, newEnergyStr, std::regex_constants::format_first_only);

    rho_norm_name = TString(fname);             // Convert back to TString
}
