
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

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
void setBin( TH1F * h,  int z_bin, double events[bins_Z][bins_p][3],  double weights, correctionTools corr, int corrType, int matchType, int chargeIdx);

int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [input k file] [input r file] [Output File]\n";
		cerr << "[Acceptance Matching Type (2,3 etc)] \n";
		cerr << "[Apply Corrections? (1 - MC, 2 - MC + pi2k, 3 - MC + pi2k + k2pi)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	double inBeam = atof(argv[1]);

	//TString in_name = argv[1];
	//TString k_name = argv[2];
    //TString r_name = argv[3];
	TString out_name = argv[2];
	int matchType = atoi(argv[3]);
	int applyCorr = atoi(argv[4]);
	int map = atoi(argv[5]);
	TString acc_name = argv[6];
	TString rho_norm_name = argv[7];

	cout<<rho_norm_name<<std::endl;
	cout<<"Create and Load files\n";

	TFile * outFile = new TFile(out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TChain * dChain = new TChain( "ePi" );
	TChain * kChain = new TChain( "ePi" );
	TChain * rChain = new TChain( "ePi" );

	//Add files
	std::cout<<"test\n";

	TString base = "/volatile/clas12/users/jphelan/SIDIS/data/";

	if( inBeam == 0 || inBeam == 10.2){
		dChain->Add( base + "final_skims/10.2/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.2/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "rho_skims/rotated_10.2_acc.root");
		}
	}
	if( inBeam == 0 || inBeam == 10.4){
		dChain->Add( base + "final_skims/10.6/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.4/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "rho_skims/rotated_10.4_acc.root");
		}
	}
	if( inBeam == 0 || inBeam == 10.6){
		dChain->Add( base + "final_skims/10.6/final_skim.root");
		if( applyCorr > 3 ){
			kChain->Add( base + "final_skims/kaons_10.6/final_skim.root");
		}
		if( applyCorr > 4 ){
			rChain->Add( base + "rho_skims/rotated_10.6_acc.root");
		}
	}

	TFile * rho_norms = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/" + rho_norm_name);
	TH3F * hNorms_pip = (TH3F*)rho_norms->Get("hNorm_pip");
	TH3F * hNorms_pim = (TH3F*)rho_norms->Get("hNorm_pim");
	TVector3 * bounds = (TVector3*)rho_norms->Get("bounds");
	double Mx_min = bounds->X();
	double Mx_max = bounds->Y();
	
	cout<<"Found rho norms\n";

	//Declare counts and weight arrays

	double events_in_bin[4][2][bins_Q2][bins_xB][bins_Z][bins_p][3] = {0}; //[sample][charge][q2][xb][z][p][Ebeam]
	double weights_in_bin[4][2][bins_Q2][bins_xB][bins_Z] = {0}; 

	TH1F * hZ[bins_Q2][bins_xB][2];
	TH1F * hZ_k[bins_Q2][bins_xB][2];
	TH1F * hZ_r[bins_Q2][bins_xB][2];
	TH1F * hZ_r_bac[bins_Q2][bins_xB][2];

	TString charge_str[2] = {"", "Pim"};

	cout<<"Make histograms\n";

	for( int i = 0; i < bins_Q2; i++ ){
		for( int j = 0; j < bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_k[i][j][k] = new TH1F("hRatio_k" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_r[i][j][k] = new TH1F("hRatio_r" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_r_bac[i][j][k] = new TH1F("hRatio_r_bac" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);

			}
		}
	}
	
	correctionTools corrector(2);
	if( matchType == 3){
		corrector.setWeightName( "corrections_10.2_3d_AN.root");
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

		//if( event_count == 1000000){break;}

		if( *eBeam != beam_energy ){
			if( matchType == 3 ) corrector.setWeightName( Form("corrections_%0.1f_3d_AN.root", *eBeam));
			else corrector.setWeightName( Form("corrections_%0.1f_AN.root", *eBeam));

			corrector.loadHistograms();	
			beam_energy = *eBeam;
		}
		TVector3 e_mom = e->get3Momentum();
		if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			double p_pi = pi[i].get3Momentum().Mag();

			if(!isGoodPion_no_acc[i]) continue;

			if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;

			int this_bin_E = (int)( (beam_energy - 10.2)/.2 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z);
			int kaon_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z);

		


			int p_bin = -1;
			for( int bin = 0; bin < 4; bin++ ){
				if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
			}

			bool matching = true;
			if( matchType == 2 ){ 
				matching = anal.applyAcceptanceMatching(pi[i], 2);
				//matching = isGoodPion[i]; }
			}
			else if( matchType == 3 ){ 
				matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), 1 ) >= 0 &&
							anal.applyAcceptanceMap( p_pi, rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), 2 ) >= 0;
				//matching = isGoodPion3d[i];
			}


			if( !matching ){ continue; }

			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );

			double weight = 1;
			double eventWeightErr = 0;

			double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
			double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
			double mc_weight = corrector.getCorrectionFactor( 2, chargeIdx );
			double k_weight = corrector.getCorrectionFactor(3, chargeIdx);
			
			
			if( applyCorr == 1 ) weight *= bin_weight;
			if( applyCorr == 2 ) weight *= acc_weight;
			if( applyCorr > 2 ){
				//weight *= acc_weight*bin_weight;
				weight *= mc_weight;
			}
			if( applyCorr > 3 ){
				if( p_bin < 0) continue;
				weight *= k_weight;
			}	
			events_in_bin[0][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin][this_bin_E]++;
			weights_in_bin[0][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]+= weight;
			
			//hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ(), weight );
			//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			

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
			if( *eBeam_k != beam_energy ){
				if( matchType == 3 ) corrector.setWeightName( Form("corrections_%0.1f_3d_AN.root", *eBeam_k));
				else corrector.setWeightName( Form("corrections_%0.1f_AN.root", *eBeam_k));
				corrector.loadHistograms();
				//corrector.loadNewEnergy( *eBeam_k );
				beam_energy = *eBeam_k;
			}
			TVector3 e_mom = e_k->get3Momentum();
			if(map && anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) <0 ) continue;
	
			for( int i = 0; i < (int) ( k.end() - k.begin() ); i++ ){
				
				int chargeIdx = (int)( k[i].getCharge() < 1 );
				double p_pi = k[i].get3Momentum().Mag();
				if(map && anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), chargeIdx + 1 ) <0)continue;


				int this_bin_E = (int)( (beam_energy - 10.2)/.2 );

				int this_bin_Q2 = (int)( ( (e_k->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
				int this_bin_xB = (int)( ( (e_k->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
				int this_bin_Z = (int)( ( (k[i].getZ() - .3)/(1.-.3) )*bins_Z);
				int kaon_bin_Z = (int)( ( (k[i].getZ() - .3)/(1.-.3) )*bins_Z);
				int p_bin = -1;
				for( int bin = 0; bin < 4; bin++ ){
					if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
				}

				bool matching = true;
				//if( matchType == 2 ){ matching = isGoodKaon[i]; }
				if( matchType == 2 ){ 
					matching = anal.applyAcceptanceMatching(k[i], 2);
					//matching = isGoodPion[i]; }
				}
				else if( matchType == 3 ){ 
					matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), 1 ) >= 0 &&
								anal.applyAcceptanceMap( p_pi, rad_to_deg*k[i].get3Momentum().Phi(), rad_to_deg*k[i].get3Momentum().Theta(), 2 ) >= 0;
					//matching = isGoodPion3d[i];
				}
				if( !matching ){ continue; }

				corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );
				double weight = 1;
				double eventWeightErr = 0;	
			
				double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
				double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
				double mc_weight = corrector.getCorrectionFactor( 2, chargeIdx );
				double k_weight = corrector.getCorrectionFactor(4, chargeIdx);
				
			
				
				if( applyCorr == 1 ) weight *= bin_weight;
				if( applyCorr == 2 ) weight *= acc_weight;
				if( applyCorr > 2 ){
				//if( !isfinite(acc_weight) || acc_err/acc_weight > .2){acc_weight = 0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
				//if( !isfinite(bin_weight) || bin_err/bin_weight > .2){bin_weight = 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
				
				//weight *= acc_weight*bin_weight;
				weight *= mc_weight;
			}
			if( applyCorr > 3 ){ 
				weight *= k_weight;
			}	
			events_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin][this_bin_E]++;
			weights_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z] += weight;

			//events_in_bin[1][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin]++;
			//hZ_k[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( k[i].getZ(), weight );
				//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
				

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
	TTreeReaderArray<bool> isGoodRho(reader_r, "isGoodPion");
	TTreeReaderValue<TLorentzVector> beam(reader_r, "beam");

	event_total = reader_r.GetEntries();
	//double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	if( applyCorr > 4 ){
		while (reader_r.Next()) {
			int event_count = reader_r.GetCurrentEntry();

			if(event_count%100000 == 0){
				cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
			}
			if( beam->E() != beam_energy ){
				if( matchType == 3 ) corrector.setWeightName( Form("corrections_%0.1f_3d_AN.root", beam->E()));
				else corrector.setWeightName( Form("corrections_%0.1f_AN.root",beam->E()));
				corrector.loadHistograms();
				beam_energy = beam->E();
			}

			TVector3 e_mom = e_r->get3Momentum();

			for( int i = 0; i < (int) ( r.end() - r.begin() ); i++ ){
				
				if( !isGoodRho[i] )continue;
				int chargeIdx = (int)( r[i].getCharge() < 1 );
				double p_pi = r[i].get3Momentum().Mag();

				int this_bin_E = (int)( (beam_energy - 10.2)/.2 );
				int this_bin_Q2 = (int)( ( (e_r->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
				int this_bin_xB = (int)( ( (e_r->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
				int this_bin_Z = (int)( ( (r[i].getZ() - .3)/(1.-.3) )*bins_Z);
				int kaon_bin_Z = (int)( ( (r[i].getZ() - .3)/(1.-.3) )*bins_Z);
				int p_bin = -1;
				for( int bin = 0; bin < 4; bin++ ){
					if( p_pi > p_bin_edges[bin] && p_pi < p_bin_edges[bin+1] ) p_bin = bin;
				}

				bool matching = true;
				//if( matchType == 2 ){ 
				//	matching = anal.applyAcceptanceMatching(pi[i], 2);
					//matching = isGoodPion[i]; }
				//}
				//if( matchType == 2 ){ matching = isGoodRho[i]; }
				//if( matchType == 2 ){ 
				//	matching = anal.applyAcceptanceMatching(r[i], 2);
					//matching = isGoodPion[i]; }
				//}
				if( matchType == 3 ){ 
					matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*r[i].get3Momentum().Phi(), rad_to_deg*r[i].get3Momentum().Theta(), 1 ) >= 0 &&
								anal.applyAcceptanceMap( p_pi, rad_to_deg*r[i].get3Momentum().Phi(), rad_to_deg*r[i].get3Momentum().Theta(), 2 ) >= 0;
					//matching = isGoodPion3d[i];
				}

				if( !matching ){ continue; }

				corrector.setKinematics( e_r->getXb(), e_r->getQ2(), r[i].getZ(), p_pi );
				double weight = 1;// rhoWeight[i];
				//if( rhoWeight[i] > 20 )continue;  
				double eventWeightErr = 0;

				double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
				double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
				double mc_weight = corrector.getCorrectionFactor( 2, chargeIdx );
				double k_weight = corrector.getCorrectionFactor( 3, chargeIdx );

			
				
				if( applyCorr == 1 ) weight *= bin_weight;
				if( applyCorr >= 2 ) weight *= acc_weight;
				if( applyCorr > 2 ){
					weight *= mc_weight;
				}
				if( applyCorr > 3 ){
					if( p_bin < 0) continue;
						weight *= k_weight;
				}	
				if( *Mx_2pi < Mx_min ){
					events_in_bin[2][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin][this_bin_E]++;
					weights_in_bin[2][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]+= weight;
					//events_in_bin[2][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin]++;
					//hZ_r[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( r[i].getZ(), weight );
				}
				if( *Mx_2pi > Mx_min && *Mx_2pi < Mx_max ){
					events_in_bin[3][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin][this_bin_E]++;
					weights_in_bin[3][chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]+= weight;
					//events_in_bin[3][chargeIdx][this_bin_Q2][this_bin_xB][kaon_bin_Z][p_bin]++;
					//hZ_r_bac[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( r[i].getZ(), weight );
				}				

			}
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
			
			for( int k = 0; k < bins_Z; k++ ){
				corrector.setKinematics( xB_min + (j-1)*.04 + 0.02, 
											Q2_min + (i-1)*.5 + 0.25,
											hZ[i-1][j-1][0]->GetXaxis()->GetBinCenter(k+1), 1.5);
				if( applyCorr <= 3){
					setBin( hZ[i-1][j-1][0],  k, events_in_bin[0][0][i-1][j-1], weights_in_bin[0][0][i-1][j-1][k], corrector, applyCorr, matchType, 0);
					setBin( hZ[i-1][j-1][1],  k, events_in_bin[0][1][i-1][j-1], weights_in_bin[0][1][i-1][j-1][k], corrector, applyCorr, matchType, 1);
				}
				else{
					setBin( hZ[i-1][j-1][0],  k, events_in_bin[0][0][i-1][j-1], weights_in_bin[0][0][i-1][j-1][k], corrector, 4, matchType, 0);
					setBin( hZ[i-1][j-1][1],  k, events_in_bin[0][1][i-1][j-1], weights_in_bin[0][1][i-1][j-1][k], corrector, 4, matchType, 1);

					setBin( hZ_k[i-1][j-1][0],  k, events_in_bin[1][0][i-1][j-1], weights_in_bin[1][0][i-1][j-1][k], corrector, 5, matchType, 0);
					setBin( hZ_k[i-1][j-1][1],  k, events_in_bin[1][1][i-1][j-1], weights_in_bin[1][1][i-1][j-1][k], corrector, 5, matchType, 1);
					
					setBin( hZ_r[i-1][j-1][0],  k, events_in_bin[2][0][i-1][j-1], weights_in_bin[2][0][i-1][j-1][k], corrector, 3, matchType, 0);
					setBin( hZ_r[i-1][j-1][1],  k, events_in_bin[2][1][i-1][j-1], weights_in_bin[2][1][i-1][j-1][k], corrector, 3, matchType, 1);

					setBin( hZ_r_bac[i-1][j-1][0],  k, events_in_bin[3][0][i-1][j-1], weights_in_bin[3][0][i-1][j-1][k], corrector, 3, matchType, 0);
					setBin( hZ_r_bac[i-1][j-1][1],  k, events_in_bin[3][1][i-1][j-1], weights_in_bin[3][1][i-1][j-1][k], corrector, 3, matchType,1);

				}

			}

			
			if( applyCorr > 3 ){
				hZ[i-1][j-1][0]->Add( hZ_k[i-1][j-1][0] );
				hZ[i-1][j-1][1]->Add( hZ_k[i-1][j-1][1] );
			}
		
			if( applyCorr > 4 ){
				for( int k = 1; k <= bins_Z; k++ ){
					double pip_cont = hZ_r_bac[i-1][j-1][0]->GetBinContent(k);
					double pip_err = hZ_r_bac[i-1][j-1][0]->GetBinError(k);
					double pip_scale = hNorms_pip->GetBinContent(j, i, k);
					double pip_scale_err = 0;//hNorms_pip->GetBinError(j, i, k);
					double pip_bin_err = pip_cont*pip_scale*sqrt( pow( pip_err/pip_cont, 2) 	+ pow( pip_scale_err/pip_scale, 2) );

					hZ_r_bac[i-1][j-1][0]->SetBinContent( k, pip_cont*pip_scale );
					if( isnan(pip_bin_err) || !isfinite(pip_bin_err)){
						cout<<"ERROR "<<std::endl;
					}
					else hZ_r_bac[i-1][j-1][0]->SetBinError( k, pip_bin_err );
					
					double pim_cont = hZ_r_bac[i-1][j-1][1]->GetBinContent(k);
					double pim_err = hZ_r_bac[i-1][j-1][1]->GetBinError(k);
					double pim_scale = hNorms_pim->GetBinContent(j, i, k);
					double pim_scale_err = 0;//hNorms_pim->GetBinError(j, i, k);
					double pim_bin_err = pim_cont*pim_scale*sqrt( pow( pim_err/pim_cont, 2) 	+ pow( pim_scale_err/pim_scale, 2) );
					
					hZ_r_bac[i-1][j-1][1]->SetBinContent( k, pim_cont*pim_scale );
					if( isnan(pim_bin_err)|| !isfinite(pim_bin_err)){
						cout<<"ERROR "<<pim_cont*pim_scale*sqrt( pow( pim_err/pim_cont, 2) 	+ pow( pim_scale_err/pim_scale, 2) )<<std::endl;
					}
					else hZ_r_bac[i-1][j-1][1]->SetBinError( k, pim_bin_err );
					

				}
			
				hZ_r[i-1][j-1][0]->Add( hZ_r_bac[i-1][j-1][0], -1 );
				hZ_r[i-1][j-1][1]->Add( hZ_r_bac[i-1][j-1][1], -1 );

				zeroSuppress(hZ_r[i-1][j-1][0]);
				zeroSuppress(hZ_r[i-1][j-1][1]);

				hZ[i-1][j-1][0]->Add( hZ_r[i-1][j-1][0], -1 );
				hZ[i-1][j-1][1]->Add( hZ_r[i-1][j-1][1], -1 );
			
			}

			hZ[i-1][j-1][0]->Divide(hZ[i-1][j-1][1]);

			TH1F * helper_3 = (TH1F *)hZ[i-1][j-1][0]->Clone();
			hZ[i-1][j-1][0]->Scale(-1.);
			hZ[i-1][j-1][0]->Add(helper_2);

			helper_3->Scale(4.);
			helper_3->Add( helper_1 );

			hZ[i-1][j-1][0]->Divide(helper_3);
			
			//hZ[i-1][j-1][0]->Print("ALL");
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

void setBin( TH1F * h,  int z_bin, double events[bins_Z][bins_p][3], double weights, correctionTools corr, int corrType, int matchType, int chargeIdx){
	
	
	double z_min = h->GetXaxis()->GetBinLowEdge(z_bin+1)+0.01;
	double z_max = h->GetXaxis()->GetBinUpEdge(z_bin+1);

	int zBinMin = (int)( ( (z_min - .3)/(1.-.3) )*bins_Z);
	int zBinMax = zBinMin+1;

	if ( weights == 0 )return;

	double error = 0;
	double term_1 = 0;
	double term_2 = 0;
	double term_3 = 0;
	double beam_energy = 10.2;

	double mc_weight = 1;//corr.getCorrectionFactor( 2, chargeIdx );
	double k_weight = 1;

	double mc_err = 0;//corr.getCorrectionError(0, chargeIdx);
	double k_err = 0;//corr.getCorrectionError(1, chargeIdx); 


	for( int E = 0; E < 3; E++){

		if( matchType == 3 ) corr.setWeightName( Form("corrections_%0.1f_3d_AN.root", (10.2 + E*.2)));
		else corr.setWeightName( Form("corrections_%0.1f_AN.root",(10.2 + E*.2)));
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
		

		double term_0 = 0;
		for( int i = 0; i < 4; i++ ){
			double p = (p_bin_edges[i] + p_bin_edges[i+1])/2.;
			for( int j = zBinMin; j < zBinMax; j++){
				
				double z = 0.3 + 0.05*j + 0.05/2.;
				corr.setKinematics(corr.getX(), corr.getQ2(), z, p);

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
			
				term_0 += k_weight*events[j][i][E];
				term_3 += pow( mc_weight*k_weight, 2)*events[j][i][E];
			}
		}
		//cout<<"TERM 0 "<<term_0<<std::endl;
		term_1 += pow( mc_err*term_0, 2);
		//cout<<"MC ERR : "<<term_1<<std::endl;
	}

	
	double term_0 = 0;
	for( int i = 0; i < 4; i++ ){
		double p = (p_bin_edges[i] + p_bin_edges[i+1])/2.;
		for( int j = zBinMin; j < zBinMax; j++){
			double z = 0.3 + 0.05*j + 0.05/2.;
			corr.setKinematics(corr.getX(), corr.getQ2(), z, p);

			switch ( corrType ){
				case 4:
					mc_weight = corr.getCorrectionFactor(2, chargeIdx);
					mc_err = corr.getCorrectionError(2, chargeIdx);
					k_weight = corr.getCorrectionFactor(3, chargeIdx);
					k_err = corr.getCorrectionError(3, chargeIdx);
					break;
				case 5:
					mc_weight = corr.getCorrectionFactor(2, chargeIdx);
					mc_err = corr.getCorrectionError(2, chargeIdx);
					k_weight = corr.getCorrectionFactor(4, chargeIdx);
					k_err = corr.getCorrectionError(4, chargeIdx);
					break;
					
			}
			double term_0 = 0;
			for( int E = 0; E < 3; E++){

				if( matchType == 3 ) corr.setWeightName( Form("corrections_%0.1f_3d_AN.root", (10.2 + E*.2)));
				else corr.setWeightName( Form("corrections_%0.1f_AN.root",(10.2 + E*.2)));
				corr.loadHistograms();
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
						k_weight = corr.getCorrectionFactor(3, chargeIdx);
						k_err = corr.getCorrectionError(3, chargeIdx);
						break;
					case 5:
						mc_weight = corr.getCorrectionFactor(2, chargeIdx);
						mc_err = corr.getCorrectionError(2, chargeIdx);
						k_weight = corr.getCorrectionFactor(4, chargeIdx);
						k_err = corr.getCorrectionError(4, chargeIdx);
						break;
				}

				term_0 += mc_weight*events[j][i][E];
				
			}
			term_2 += pow( term_0*k_err, 2);
			

		}
	}
	
	error = sqrt( term_1 + term_2 + term_3);

	h->SetBinContent( z_bin +1, weights);
	h->SetBinError(z_bin + 1, error);
}