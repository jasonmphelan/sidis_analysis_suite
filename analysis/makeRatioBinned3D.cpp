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
TString var_name = "null";
int bins_var = 1;
double var_min = 0;
double var_max = 1;

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
		cerr << "       [match type] [apply corr] [map] [acceptance file]\n";
		cerr << "       [var name] [bins var] [var min] [var max] [binned out] [correction name]\n";
		cerr << "       (opt: Mx2 cut) (opt: rho mode)\n";
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

	// Probe correction file for 4D binning metadata before any array allocation
	if( bins_var > 1 ){
		TString probePath = (TString)_DATA + "/correctionFiles/" + correction_name;
		TFile * pf = TFile::Open(probePath);
		if( pf && !pf->IsZombie() ){
			int count = 0;
			while( pf->Get(Form("hMcCorrP_%d", count)) ) count++;
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
	

	if( inBeam == 0 || inBeam == 10.2){
		dChain->Add( base + "final_skims/10.2/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.2/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "final_skims/rho_skims/rotated_10.2_sym.root");
		}
	}
	if( inBeam == 0 || inBeam == 10.4){
		dChain->Add( base + "final_skims/10.4/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.4/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "final_skims/rho_skims/rotated_10.4_sym.root");
			
		}
	}
	if( inBeam == 0 || inBeam == 10.6){
		dChain->Add( base + "final_skims/10.6/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.6/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "final_skims/rho_skims/rotated_10.6_sym.root");

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
	corrector.setWeightName( correction_name );
	corrector.setWeightName4D( correction_name );
	if( n4d_corr_bins > 0 ){
		corrector.setN4dBins(n4d_corr_bins);
	}
	corrector.setK2piName( "corrections_k2pi_AN.root");
	corrector.setPi2kName( "corrections_pi2k_AN.root");
	corrector.loadHistograms();	
	corrector.printFilePaths();

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions("matchCut2D_map.root");
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

	// Lambda: compute all bin indices for a given electron/pion/momentum
	struct Bins { int E, Q2, xB, Z, var, p; };
	auto calcBins = [&](const electron& el, const pion& pi_ref, double p_mom) -> Bins {
		int bE   = (int)((beam_energy - 10.2) / .2);
		int bQ2  = (int)(((el.getQ2()  - Q2_min) / (Q2_max - Q2_min)) * bins_Q2);
		int bxB  = (int)(((el.getXb()  - xB_min) / (xB_max - xB_min)) * bins_xB);
		int bZ   = (int)(((pi_ref.getZ() - .3)    / (1. - .3))         * bins_Z);
		int bVar = (int)(((getVarVal(var_name, el, pi_ref) - var_min) / (var_max - var_min)) * bins_var);
		int bP   = -1;
		for (int bin = 0; bin < 4; bin++)
			if (p_mom > p_bin_edges[bin] && p_mom < p_bin_edges[bin+1]) { bP = bin; break; }
		return {bE, bQ2, bxB, bZ, bVar, bP};
	};

	// Lambda: apply acceptance matching based on matchType
	auto applyMatching = [&](const pion& pi_ref, double p_mom) -> bool {
		if (matchType == 1 || matchType == 2)
			return anal.applyAcceptanceMatching(pi_ref, 2);
		if (matchType == 3)
			return anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 1) >= 0 &&
			       anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 2) >= 0;
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
		if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
		
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			double p_pi = pi[i].get3Momentum().Mag();

			if(!isGoodPion_no_acc[i]) continue;
			if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

			if( !applyMatching(pi[i], p_pi) ) continue;

			auto b = calcBins(*e, pi[i], p_pi);
			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );
			double weight = 1;
			if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
			if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( applyCorr > 2 )  weight *= corrector.getCorrectionFactor(2, chargeIdx);
			if( applyCorr > 3 && b.p >= 0) weight *= corrector.getCorrectionFactor(3, chargeIdx);
			events_in_bin[0][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.p][b.E]++;
			weights_in_bin[0][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.E] += weight;
			errors_in_bin[0][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.p][b.E]++;

		}
	}
	
	///////////////////////////////////////
	////////////// Kaons /////////////////
	///////////////////////////////////////

	//TTreeReader reader_k("ePi", kFile);
	TTreeReader reader_k(kChain);
	TTreeReaderValue<double> eBeam_k( reader_k, "Ebeam" );
	TTreeReaderValue<electron> e_k(reader_k, "e");
	TTreeReaderArray<pion> k(reader_k, "pi");
	TTreeReaderArray<bool> isGoodKaon_no_acc(reader_k, "isGoodPion_no_acc");


	event_total = reader_k.GetEntries();

	if( applyCorr > 3 ){
		cout << "Begin kaon analysis\n";
		while (reader_k.Next()) {
			int event_count = reader_k.GetCurrentEntry();
			if(event_count%100000 == 0){
				cout << "Kaons: " << event_count << " / " << event_total << "\n";
			}
			
			updateCorrectionsForBeam(*eBeam_k, beam_energy, correction_name, corrector);

			TVector3 e_mom = e_k->get3Momentum();
			if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
	
			for( int i = 0; i < (int) ( k.end() - k.begin() ); i++ ){
				
				int chargeIdx = (int)( k[i].getCharge() < 1 );
				double p_pi = k[i].get3Momentum().Mag();
				if( !isGoodKaon_no_acc[i] )continue;

				if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

				if( !applyMatching(k[i], p_pi) ) continue;

				auto b = calcBins(*e_k, k[i], p_pi);
				corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );
				if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
				double weight = 1;
				if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
				if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
				if( applyCorr > 2 )  weight *= corrector.getCorrectionFactor(2, chargeIdx);
				if( applyCorr > 3 )  weight *= corrector.getCorrectionFactor(4, chargeIdx);
				events_in_bin[1][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.p][b.E]++;
				weights_in_bin[1][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.E] += weight;
				errors_in_bin[1][chargeIdx][b.var][b.Q2][b.xB][b.Z][b.p][b.E]++;

			}
		}
	}
	
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

	if( applyCorr > 4 ){
		cout << "Begin rho analysis\n";
		while (reader_r.Next()) {
			int event_count = reader_r.GetCurrentEntry();

			if(event_count%100000 == 0){
				cout << "Rhos:  " << event_count << " / " << event_total << "\n";
			}
			updateCorrectionsForBeam(beam->E(), beam_energy, correction_name, corrector);
			
			
			TVector3 e_mom = e_r->get3Momentum();

			int i = -1; 
			

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
			corrector.setKinematics( e_r->getXb(), e_r->getQ2(), r[i].getZ(), p_pi );
			if( n4d_corr_bins > 0 ) corrector.set4dBin(b.var);
			double weight = rhoWeight[i];
			if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( applyCorr > 2 )  weight *= corrector.getCorrectionFactor(2, chargeIdx);
			if( applyCorr > 3 && b.p >= 0) weight *= corrector.getCorrectionFactor(3, chargeIdx);

			events_in_bin[index][0][b.var][b.Q2][b.xB][b.Z][b.p][b.E] += rhoWeight[i];
			weights_in_bin[index][0][b.var][b.Q2][b.xB][b.Z][b.E] += weight;
			errors_in_bin[index][0][b.var][b.Q2][b.xB][b.Z][b.p][b.E] += pow(rhoWeight[i],2) + pow(rhoError[i],2);
			events_in_bin[index][1][b.var][b.Q2][b.xB][b.Z][b.p][b.E] += rhoWeight[i];
			weights_in_bin[index][1][b.var][b.Q2][b.xB][b.Z][b.E] += weight;
			errors_in_bin[index][1][b.var][b.Q2][b.xB][b.Z][b.p][b.E] += pow(rhoWeight[i],2) + pow(rhoError[i],2);
		}
	}




	//////////////////////////////////////
	///////////// Calcs //////////////////
	//////////////////////////////////////

	cout<<"Start calculations\n";



	for( int i = 1; i <= bins_Q2; i++ ){
		for( int j = 1; j <= bins_xB; j++ ){
			for( int var = 0; var < bins_var; var++ ){
				for( int E = 0; E < 3; E++ ){	
					
					for( int k = 0; k < bins_Z; k++ ){
						
						corrector.setKinematics( xB_min + (j-1)*.04 + 0.02, 
													Q2_min + (i-1)*.5 + 0.25,
													hZ[var][i-1][j-1][0]->GetXaxis()->GetBinCenter(k+1), 1.5);
						if( n4d_corr_bins > 0 ) corrector.set4dBin(var);

						if( applyCorr <= 3){
							setBin( hZ[var][i-1][j-1][0],  k, events_in_bin[0][0][var][i-1][j-1], errors_in_bin[0][0][var][i-1][j-1], weights_in_bin[0][0][var][i-1][j-1][k][E], corrector, applyCorr,  0, E, correction_name );
							setBin( hZ[var][i-1][j-1][1],  k, events_in_bin[0][1][var][i-1][j-1], errors_in_bin[0][1][var][i-1][j-1], weights_in_bin[0][1][var][i-1][j-1][k][E], corrector, applyCorr, 1, E, correction_name);
						}
						else{
							setBin( hZ[var][i-1][j-1][0],  k, events_in_bin[0][0][var][i-1][j-1], errors_in_bin[0][0][var][i-1][j-1], weights_in_bin[0][0][var][i-1][j-1][k][E], corrector, 4, 0, E, correction_name );
							setBin( hZ[var][i-1][j-1][1],  k, events_in_bin[0][1][var][i-1][j-1], errors_in_bin[0][1][var][i-1][j-1],weights_in_bin[0][1][var][i-1][j-1][k][E], corrector, 4, 1, E, correction_name );

							setBin( hZ_k[var][i-1][j-1][0],  k, events_in_bin[1][0][var][i-1][j-1], errors_in_bin[1][0][var][i-1][j-1],weights_in_bin[1][0][var][i-1][j-1][k][E], corrector, 5, 0, E, correction_name );
							setBin( hZ_k[var][i-1][j-1][1],  k, events_in_bin[1][1][var][i-1][j-1], errors_in_bin[1][1][var][i-1][j-1],weights_in_bin[1][1][var][i-1][j-1][k][E], corrector, 5, 1, E, correction_name );
							
							setBin( hZ_r[var][i-1][j-1][0],  k, events_in_bin[2][0][var][i-1][j-1], errors_in_bin[2][0][var][i-1][j-1],weights_in_bin[2][0][var][i-1][j-1][k][E], corrector, 4, 0, E, correction_name );
							setBin( hZ_r[var][i-1][j-1][1],  k, events_in_bin[2][1][var][i-1][j-1], errors_in_bin[2][1][var][i-1][j-1],weights_in_bin[2][1][var][i-1][j-1][k][E], corrector, 4, 1, E, correction_name );

							
							

						}

					}

				}

				outFile->cd();
				hZ[var][i-1][j-1][0]->Write();
				hZ[var][i-1][j-1][1]->Write();
				if( applyCorr > 3 ){
					hZ_k[var][i-1][j-1][0]->Write();
					hZ_k[var][i-1][j-1][1]->Write();
				}
				if( applyCorr > 4 ){
					hZ_r[var][i-1][j-1][0]->Write();
					hZ_r[var][i-1][j-1][1]->Write();
				}
			}
		}
	}

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

    corrector.loadHistograms();
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
    

    corrector.loadHistograms();
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
