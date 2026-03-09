
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <regex>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TChain.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "TRandom3.h"
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

void zeroSuppress( TH1F * h);
double getVarVal( electron e, pion pi ); 
void setBin( TH1F * h,  int z_bin, double events[bins_Z][bins_p][3],  double errors[bins_Z][bins_p][3], double weights, correctionTools corr, int corrType, int matchType, int chargeIdx, int E_bin, double norm = 1);
void ScaleHistogram(TH1* h, const TH3* hNorms, int j, int i, bool includeScaleError = false);
bool updateCorrectionsForBeam(
    double eBeam,
    double& beam_energy,
    int matchType,
    correctionTools& corrector
);

void getNewRhoName( double E, TString& rho_norm_name) ;



int main(int argc, char** argv) {

    // Minimum required arguments:
    // argv[0] = program name
    // argv[1] = beam energy
    // argv[2] = output file
    // argv[3] = match type
    // argv[4] = apply corrections
    // argv[5] = map
    // argv[6] = acceptance file
    // argv[7] = rho file
    // argv[8] = rho norm file
    if (argc < 8) {
        cerr << "Incorrect number of arguments. Please use:\n";
        cerr << "./code [Beam Energy]\n";
        cerr << "      [Output File]\n";
        cerr << "      [Acceptance Matching Type]\n";
        cerr << "      [Apply Corrections]\n";
        cerr << "      [Map]\n";
        cerr << "      [Acceptance File]\n";
        cerr << "      [Rho File]\n";
        cerr << "      [Rho Norm File]\n";
        cerr << "      (Optional) [Extra TString 1]\n";
        cerr << "      (Optional) [Extra TString 2]\n";
        return -1;
    }

    // Required arguments
    double inBeam     = atof(argv[1]);
    TString out_name  = argv[2];
    int matchType     = atoi(argv[3]);
    int applyCorr     = atoi(argv[4]);
    int map           = atoi(argv[5]);
    TString acc_name  = argv[6];
	double Mx_cut = atof(argv[7]);

    // Optional arguments
	TString mc_file = (argc > 8)  ? TString(argv[8])  : "corrections_10.2_AN_test.root";
    TString pi2k_file = (argc > 9)  ? TString(argv[9])  : "corrections_pi2k_AN.root";
    TString k2pi_file = (argc > 10) ? TString(argv[10]) : "corrections_k2pi_AN.root";

    cerr << "Files used:\n";
    cerr << "  Beam Energy     : " << inBeam << "\n";
    cerr << "  Output File     : " << out_name << "\n";
    cerr << "  Acc File        : " << acc_name << "\n";
    cerr << "  Rho File        : " << Mx_cut << "\n";
	
	if (!mc_file.IsNull())
		cerr << "  Extra Arg 0     : " <<  mc_file << "\n";
    if (!pi2k_file.IsNull())
        cerr << "  Extra Arg 1     : " <<  pi2k_file << "\n";
    if (!k2pi_file.IsNull())
        cerr << "  Extra Arg 2     : " << k2pi_file << "\n";

    // --- rest of your code ---


	TFile * outFile = new TFile(out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TChain * dChain = new TChain( "ePi" );
	TChain * kChain = new TChain( "ePi" );
	TChain * rChain = new TChain( "ePi" );

	//Add files
	TString base = "../trees/"; //"/volatile/clas12/users/jphelan/SIDIS/data/";
	constexpr int N_E = 3;

	
	double     E_list[N_E]       = {10.2, 10.4, 10.6};  // or fill later


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



	cout<<"Found rho norms\n";

	//Declare counts and weight arrays

	double events_in_bin[4][2][bins_Q2][bins_xB][bins_Z][bins_p][3] = {0}; //[sample][charge][q2][xb][z][p][Ebeam]
	double weights_in_bin[4][2][bins_Q2][bins_xB][bins_Z][3] = {0}; 
	double errors_in_bin[4][2][bins_Q2][bins_xB][bins_Z][bins_p][3] = {0};

	TH1F * hZ[bins_Q2][bins_xB][2];
	TH1F * hZ_k[bins_Q2][bins_xB][2];
	TH1F * hZ_r[bins_Q2][bins_xB][2];

	TString charge_str[2] = {"", "Pim"};

	cout<<"Make histograms\n";

	for( int i = 0; i < bins_Q2; i++ ){
		for( int j = 0; j < bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_k[i][j][k] = new TH1F("hRatio_k" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_r[i][j][k] = new TH1F("hRatio_r" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);

			}
		}
	}
	
	int corrector_mode = (matchType==1) ? 4 : 2;
	correctionTools corrector(corrector_mode);
	if( matchType == 3){
		corrector.setWeightName( "corrections_10.2_3d_AN.root");
	}
	if( matchType == 0){
		corrector.setWeightName( "corrections_10.2_no_match_AN.root");
	}
	corrector.setK2piName( k2pi_file );
	corrector.setPi2kName( pi2k_file );
	if( matchType == 1){
		corrector.setPi2kName( "corrections_10.2_no_match_AN_4D.root");
	}
	corrector.setWeightName( "corrections_10.2_AN_test.root");
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

	cout<<"Begin analysis\n";
	double beam_energy = 10.2; //current energy of file

	//TTreeReader reader_rec("ePi", inFile);
	TTreeReader reader_rec(dChain);
	TTreeReaderValue<double> eBeam( reader_rec, "Ebeam" );
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");

	TTreeReaderArray<bool> isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");

	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");

	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();
		
		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
		}

		updateCorrectionsForBeam(*eBeam, beam_energy, matchType, corrector);

		TVector3 e_mom = e->get3Momentum();
		if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
		
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			double p_pi = pi[i].get3Momentum().Mag();

			if(!isGoodPion_no_acc[i]) continue;
			//if( pi[i].getMx() < 1.6 )continue;
			if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

			int this_bin_E = (int)( (beam_energy - 10.2)/.2 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z);

			int p_bin = -1;
			for( int bin = 0; bin < 4; bin++ ){
				if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
			}

			bool matching = true;
			double weight = 1;

			if( matchType == 1 ){
				int pt_bin = (int)( ( (pi[i].getPi_q().Pt() - 0)/(1.2-0) )*4);
				if ( pt_bin > 3 ){continue;}
				corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), pi[i].getPi_q().Pt() );
				weight *= corrector.getCorrectionFactor(3, chargeIdx);

				
			}
			if( matchType == 2 ){ matching = anal.applyAcceptanceMatching(pi[i], 2);}
			else if( matchType == 3 ){ 
				matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), 1 ) >= 0 &&
							anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), 2 ) >= 0;
			}

			if( !matching ){ continue; }

			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );

			if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( applyCorr > 2 ) weight *= corrector.getCorrectionFactor( 2, chargeIdx );
			if( applyCorr > 3  && p_bin >= 0) weight *= corrector.getCorrectionFactor(3, chargeIdx);
			events_in_bin[0][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]++;
			weights_in_bin[0][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_E]+= weight;
			errors_in_bin[0][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]++;			

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
	TTreeReaderArray<bool> isGoodKaon(reader_k, "isGoodPion");
	TTreeReaderArray<bool> isGoodKaon3d(reader_k, "isGoodPion_3d");

	event_total = reader_k.GetEntries();
	//double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	if( applyCorr > 3 ){
		while (reader_k.Next()) {
			int event_count = reader_k.GetCurrentEntry();

			if(event_count%100000 == 0){
				cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
			}
			
			updateCorrectionsForBeam(*eBeam_k, beam_energy, matchType, corrector);

			TVector3 e_mom = e_k->get3Momentum();
			if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
	
			for( int i = 0; i < (int) ( k.end() - k.begin() ); i++ ){
				
				int chargeIdx = (int)( k[i].getCharge() < 1 );
				double p_pi = k[i].get3Momentum().Mag();
				if( !isGoodKaon_no_acc[i] )continue;
				if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

				int this_bin_E = (int)( (beam_energy - 10.2)/.2 );

				int this_bin_Q2 = (int)( ( (e_k->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
				int this_bin_xB = (int)( ( (e_k->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
				int this_bin_Z = (int)( ( (k[i].getZ() - .3)/(1.-.3) )*bins_Z);
				int p_bin = -1;
				for( int bin = 0; bin < 4; bin++ ){
					if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
				}

				bool matching = true;
				double weight = 1;
				

				if( matchType == 1 ){
					int pt_bin = (int)( ( (k[i].getPi_q().Pt() - 0)/(1.2-0) )*4);
					if ( pt_bin > 3 ){continue;}
					corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );
					weight *= corrector.getCorrectionFactor(3, chargeIdx);
				}
				
				if( matchType == 2 ){ 
					matching = anal.applyAcceptanceMatching(k[i], 2);
			
				}
				else if( matchType == 3 ){ 
					matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), 1 ) >= 0 &&
								anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), 2 ) >= 0;
					
				}
				if( !matching ){ continue; }

				corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );
				
				if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
				if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
				if( applyCorr > 2 ) weight *= corrector.getCorrectionFactor( 2, chargeIdx );
			
			if( applyCorr > 3 ){ 
				weight *= corrector.getCorrectionFactor(4, chargeIdx);
			}	
			events_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]++;
			weights_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_E] += weight;
			errors_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E] ++;
				

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
	TTreeReaderArray<double> rhoWeight_sym( reader_r, "rhoWeight_sym");
	TTreeReaderValue<double> Mx_2pi( reader_r, "Mx_2pi");
	TTreeReaderArray<double> rhoError( reader_r, "rhoErr");
	//TTreeReaderArray<double> rhoError_sym( reader_r, "rhoErr_sym");
	TTreeReaderArray<bool> isGoodRho(reader_r, "isGoodPion");
	TTreeReaderValue<TLorentzVector> beam(reader_r, "beam");

	event_total = reader_r.GetEntries();

	if( applyCorr > 4 ){
		while (reader_r.Next()) {
			int event_count = reader_r.GetCurrentEntry();

			if(event_count%100000 == 0){
				cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
			}
			int this_bin_E = (int)( (beam_energy - 10.2)/.2 );
			updateCorrectionsForBeam(beam->E(), beam_energy, matchType, corrector);
			
			
			TVector3 e_mom = e_r->get3Momentum();

			int i = -1; 
			
			if( isGoodRho[1] ) i = 1; //By default, use pim
			else if( !isGoodRho[1] &&  !anal.applyAcceptanceMatching(r[1], 2) && isGoodRho[0]) i=0;
			else { continue; }

			int chargeIdx = (int)( r[i].getCharge() < 1 );
			double p_pi = r[i].get3Momentum().Mag();

			
			int this_bin_Q2 = (int)( ( (e_r->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e_r->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (r[i].getZ() - .3)/(1.-.3) )*bins_Z);
			
			int p_bin = -1;
			for( int bin = 0; bin < 4; bin++ ){
				if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
			}

			bool matching = true;
			
			if( matchType == 2 ){ 
				matching = anal.applyAcceptanceMatching(r[i], 2);
			}
			if( matchType == 3 ){ 
				matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*r[i].get3Momentum().Phi(), rad_to_deg*r[i].get3Momentum().Theta(), 1 ) >= 0 &&
							anal.applyAcceptanceMap( p_pi, rad_to_deg*r[i].get3Momentum().Phi(), rad_to_deg*r[i].get3Momentum().Theta(), 2 ) >= 0;
			}

			if( !matching ){ continue; }

			corrector.setKinematics( e_r->getXb(), e_r->getQ2(), r[i].getZ(), p_pi );
			
			//if( rhoWeight[i] > 100 || rhoWeight[i] < 0)continue;  
			double weight = rhoWeight[i];//0.5*rhoWeight[i];// - (double)(i==1);
			
		

			if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
			if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( applyCorr > 2 ){
				weight *=  corrector.getCorrectionFactor( 2, chargeIdx );
			}
			if( applyCorr > 3 && p_bin >= 0) weight *= corrector.getCorrectionFactor( 3, chargeIdx );


			if( (*Mx_2pi)*(*Mx_2pi) > Mx_cut ) continue;
			//int index = ( *Mx_2pi > Mx_min && *Mx_2pi < Mx_max) ? 3 : 2;
			int index = 2;

			events_in_bin[index][0][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+=rhoWeight[i];
			weights_in_bin[index][0][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_E]+=  weight;
			errors_in_bin[index][0][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+= pow(rhoWeight[i],2) + pow(rhoError[i],2);
						
			events_in_bin[index][1][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+=rhoWeight[i];
			weights_in_bin[index][1][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_E]+=  weight;
			errors_in_bin[index][1][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+= pow(rhoWeight[i],2) + pow(rhoError[i],2);
							

/*
				
				if( rhoWeight_sym[i] > 100 || rhoWeight_sym[i] < 0)continue;  
				weight = 0.5*rhoWeight_sym[i];// - (double)(i==1);
				double rhoError_sym[2] = {0.5*rhoError_sym[0], 0.5*rhoError_sym[1]};
				
				if( applyCorr == 1 ) weight *= corrector.getCorrectionFactor(0, (int)(!chargeIdx));
				if( applyCorr == 2 ) weight *= corrector.getCorrectionFactor(1, (int)(!chargeIdx));
				if( applyCorr > 2 ){
					weight *=  corrector.getCorrectionFactor( 2, (int)(!chargeIdx) );
				}
				if( applyCorr > 3 && p_bin >= 0) weight *= corrector.getCorrectionFactor( 3, (int)(!chargeIdx) );

				events_in_bin[index][(int)(!chargeIdx)][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+=0.5*rhoWeight_sym[i];
				weights_in_bin[index][(int)(!chargeIdx)][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_E]+= weight;
				errors_in_bin[index][(int)(!chargeIdx)][this_bin_Q2][this_bin_xB][this_bin_Z][p_bin][this_bin_E]+= pow(0.5*rhoWeight_sym[i],2) + pow(0.5*rhoError_sym[i],2);
			
*/
			
		}
	}




	//////////////////////////////////////
	///////////// Calcs //////////////////
	//////////////////////////////////////

	cout<<"Start calculations\n";

	TH1F * helper_1 = new TH1F("helper1", "helper1", 14, .3, 1);
	TH1F * helper_2 = new TH1F("helper2", "helper2", 14, .3, 1);
		

	
	for( int i = 1; i <= bins_Z; i++ ){
		helper_1->SetBinContent(i, -1.);
		helper_2->SetBinContent(i, 4.);
		
		helper_1->SetBinError(i, 0.);
		helper_2->SetBinError(i, 0.);
	}


	for( int i = 1; i <= bins_Q2; i++ ){
		for( int j = 1; j <= bins_xB; j++ ){
			for( int E = 0; E < 3; E++ ){	
				
				for( int k = 0; k < bins_Z; k++ ){
					
					corrector.setKinematics( xB_min + (j-1)*.04 + 0.02, 
												Q2_min + (i-1)*.5 + 0.25,
												hZ[i-1][j-1][0]->GetXaxis()->GetBinCenter(k+1), 1.5);
					if( applyCorr <= 3){
						setBin( hZ[i-1][j-1][0],  k, events_in_bin[0][0][i-1][j-1], errors_in_bin[0][0][i-1][j-1], weights_in_bin[0][0][i-1][j-1][k][E], corrector, applyCorr, matchType, 0, E );
						setBin( hZ[i-1][j-1][1],  k, events_in_bin[0][1][i-1][j-1], errors_in_bin[0][1][i-1][j-1], weights_in_bin[0][1][i-1][j-1][k][E], corrector, applyCorr, matchType, 1, E );
					}
					else{
						setBin( hZ[i-1][j-1][0],  k, events_in_bin[0][0][i-1][j-1], errors_in_bin[0][0][i-1][j-1], weights_in_bin[0][0][i-1][j-1][k][E], corrector, 4, matchType, 0, E );
						setBin( hZ[i-1][j-1][1],  k, events_in_bin[0][1][i-1][j-1], errors_in_bin[0][1][i-1][j-1],weights_in_bin[0][1][i-1][j-1][k][E], corrector, 4, matchType, 1, E );

						setBin( hZ_k[i-1][j-1][0],  k, events_in_bin[1][0][i-1][j-1], errors_in_bin[1][0][i-1][j-1],weights_in_bin[1][0][i-1][j-1][k][E], corrector, 5, matchType, 0, E );
						setBin( hZ_k[i-1][j-1][1],  k, events_in_bin[1][1][i-1][j-1], errors_in_bin[1][1][i-1][j-1],weights_in_bin[1][1][i-1][j-1][k][E], corrector, 5, matchType, 1, E );
						
						setBin( hZ_r[i-1][j-1][0],  k, events_in_bin[2][0][i-1][j-1], errors_in_bin[2][0][i-1][j-1],weights_in_bin[2][0][i-1][j-1][k][E], corrector, 4, matchType, 0, E );
						setBin( hZ_r[i-1][j-1][1],  k, events_in_bin[2][1][i-1][j-1], errors_in_bin[2][1][i-1][j-1],weights_in_bin[2][1][i-1][j-1][k][E], corrector, 4, matchType, 1, E );

						
					}

				}

			}
			
			outFile->cd();

			
			//hZ_k[i-1][j-1][0]->Write();
			//hZ_k[i-1][j-1][1]->Write();
			
			//hZ_r[i-1][j-1][0]->Write();
			//hZ_r[i-1][j-1][1]->Write();

			//hZ_r_bac[i-1][j-1][0]->Write();
			//hZ_r_bac[i-1][j-1][1]->Write();
			
			if( applyCorr > 3 ){
				hZ[i-1][j-1][0]->Add( hZ_k[i-1][j-1][0] );
				hZ[i-1][j-1][1]->Add( hZ_k[i-1][j-1][1] );
			}
			//hZ[i-1][j-1][0]->Write();
			//hZ[i-1][j-1][1]->Write();
			if( applyCorr > 4 ){
				
				//hZ_r[i-1][j-1][0]->Add( hZ_r_bac[i-1][j-1][0], -1 );
				//hZ_r[i-1][j-1][1]->Add( hZ_r_bac[i-1][j-1][1], -1 );

				//zeroSuppress(hZ_r[i-1][j-1][0]);
				//zeroSuppress(hZ_r[i-1][j-1][1]);

				hZ[i-1][j-1][0]->Add( hZ_r[i-1][j-1][0], -1 );
				hZ[i-1][j-1][1]->Add( hZ_r[i-1][j-1][1], -1 );
			
			}
			//hZ_r[i-1][j-1][0]->Write();
			//hZ_r[i-1][j-1][1]->Write();


			hZ[i-1][j-1][0]->Divide(hZ[i-1][j-1][1]);
			
			TH1F * helper_3 = (TH1F *)hZ[i-1][j-1][0]->Clone();
			hZ[i-1][j-1][0]->Scale(-1.);
			hZ[i-1][j-1][0]->Add(helper_2);

			helper_3->Scale(4.);
			helper_3->Add( helper_1 );

			hZ[i-1][j-1][0]->Divide(helper_3);
			
			hZ[i-1][j-1][0]->SetDirectory(outFile);
			hZ[i-1][j-1][1]->SetDirectory(outFile);

			outFile->cd();
			hZ[i-1][j-1][0]->Write();
			hZ[i-1][j-1][1]->Write();
		
		}
	}		

	cout<<"Closing Out File\n";
	delete outFile;
	cout<<"Done \n";

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


double getVarVal(TString var,  electron * e, pion pi ){
	if( var == "p_e" ) return e->get3Momentum().Mag();
	if( var == "theta_e" ) return e->get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e->get3Momentum().Phi()*rad_to_deg;
	if( var == "sector_e" ) { double phi_deg = e->get3Momentum().Phi()*rad_to_deg + 20.; if (phi_deg < 0.) phi_deg += 360.; return (int)(phi_deg / 60.) + 1; }
	if( var == "W2" ) return e->getW2();
	if( var == "Q2" ) return e->getQ2();
	if( var == "xB" ) return e->getXb();
	if( var == "y" ) return e->getY();


	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) return pi.getPi_q().Phi()*rad_to_deg;
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	return 0;
}

void setBin( TH1F * h,  int z_bin, double events[bins_Z][bins_p][3], double errors[bins_Z][bins_p][3], double weights, correctionTools corr, int corrType, int matchType, int chargeIdx, int E_bin, double norm){
	
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

	
	if( matchType == 3 ) corr.setWeightName( Form("corrections_%0.1f_3d_AN.root", (10.2 + E*.2)));
	else if( matchType == 0) corr.setWeightName( Form("corrections_%0.1f_no_match_AN.root",(10.2 + E*.2)));
	else if( matchType == 1) corr.setPi2kName( Form("corrections_%0.1f_no_match_AN_4D.root",(10.2 + E*.2)));
	else corr.setWeightName( Form("corrections_%0.1f_AN_test.root", (10.2 + E*.2)));

	corr.loadHistograms();
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
		case 3:
			mc_weight = corr.getCorrectionFactor(2, chargeIdx);
			mc_err = corr.getCorrectionError(2, chargeIdx);
			break;
		case 4:
			mc_weight = corr.getCorrectionFactor(2, chargeIdx);
			mc_err = corr.getCorrectionError(2, chargeIdx);
			break;
		case 5:
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

	h->SetBinContent( z_bin +1, bin_old + norm*weights);
	h->SetBinError(z_bin + 1, sqrt( pow(err_old, 2)+ pow(norm*error,2) ));
}


void ScaleHistogram(
    TH1* h,
    const TH3* hNorms,   // or TH2/THn depending on your object
    int j,
    int i,
    bool includeScaleError 
){
    if (!h || !hNorms) return;

    const int nBins = h->GetNbinsX();

    for (int k = 1; k <= nBins; ++k) {

        const double cont = h->GetBinContent(k);
        const double err  = h->GetBinError(k);

        const double scale     = hNorms->GetBinContent(j, i, k);
        const double scale_err = includeScaleError
                                 ? hNorms->GetBinError(j, i, k)
                                 : 0.0;

        // Set scaled content
        h->SetBinContent(k, cont * scale);

        // Protect against invalid divisions
        if (cont <= 0.0 || scale <= 0.0) continue;

        const double rel_err2 =
            std::pow(err / cont, 2) +
            std::pow(scale_err / scale, 2);

        const double bin_err = cont * scale * std::sqrt(rel_err2);

        if (std::isfinite(bin_err))
            h->SetBinError(k, bin_err);
    }
}


bool updateCorrectionsForBeam(
    double eBeam,
    double& beam_energy,
    int matchType,
    correctionTools& corrector
){
    if (eBeam == beam_energy)
        return false;

    corrector.setWeightName(
        Form(matchType == 3
             ? "corrections_%0.1f_3d_AN.root"
             : "corrections_%0.1f_AN_test.root",
             eBeam)
    );
	if( matchType == 0){
		corrector.setWeightName(Form("corrections_%0.1f_no_match_AN.root",eBeam));
	}

	

    corrector.loadHistograms();
    beam_energy = eBeam;
    return true;
}


// Helper function to load RhoNorms for multiple beam energies


//void getNewRhoName(double E, TString& rho_norm_name) {
//	int fileNum = 10*(E-10);/
//	rho_norm_name.ReplaceAll("2", E);
//	rho_norm_name.ReplaceAll("4", E);
//	rho_norm_name.ReplaceAll("6", E);
//}
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
