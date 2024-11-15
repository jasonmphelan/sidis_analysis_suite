
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
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST


using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;


int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [input k file] [Output File]\n";
		cerr << "[Acceptance Matching Type (2,3 etc)] \n";
		cerr << "[Acc/Bin Correction?]\n";
		cerr << "[Kaon Correction?]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
	TString k_name = argv[2];
       	TString out_name = argv[3];
	int matchType = atoi(argv[4]);
	int applyCorr = atoi(argv[5]);
	int applyKaon;
       	
	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TFile * inFile = new TFile( in_name );
	TFile * kFile = new TFile( k_name );

	TH1F * hZ_uw[bins_Q2+1][bins_xB+1][2];
	TH1F * hZ_k_uw[bins_Q2+1][bins_xB+1][2];
	
	TH1F * hZ[bins_Q2+1][bins_xB+1][2];
	TH1F * hZ_k[bins_Q2+1][bins_xB+1][2];

	TString charge_str[2] = {"", "Pim"};

	for( int i = 0; i < bins_Q2; i++ ){
		for( int j = 0; j < bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
				hZ_uw[i][j][k] = new TH1F("hRatio_uw" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_k_uw[i][j][k] = new TH1F("hRatio_k_uw" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
				hZ_k[i][j][k] = new TH1F("hRatio_k" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
						
			}
		}
	}
	
	correctionTools corrector(1);
	corrector.loadHistograms();	

	TTreeReader reader_rec("ePi", inFile);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");

	int event_total = reader_rec.GetEntries();
	double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
		}

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			double p_pi = pi[i].get3Momentum().Mag();
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z);

			bool matching = true;
			if( matchType == 2 ){ matching = !isGoodPion[i]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[i]; }
			//else{ matching = false; }

			if( matching ){ continue; }
			/*
			double acc_weight = 1;
			double bin_weight = 1;
			double acc_err = 0;
			double bin_err = 0;	
			
			if( applyCorr == 1 ){
				acc_weight = (double) accWeight[chargeIdx]->GetBinContent( this_bin_xB, this_bin_Q2, this_bin_Z );
				bin_weight = (double) binWeight[chargeIdx]->GetBinContent( this_bin_xB, this_bin_Q2, this_bin_Z );
				acc_err = (double) accWeight[chargeIdx]->GetBinError( this_bin_xB, this_bin_Q2, this_bin_Z );
				bin_err = (double) binWeight[chargeIdx]->GetBinError( this_bin_xB, this_bin_Q2, this_bin_Z );
				
				if( !isfinite(acc_weight) || acc_err/acc_weight > .2){acc_weight = 0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
				if( !isfinite(bin_weight) || bin_err/bin_weight > .2){bin_weight = 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
			}
			*/
				
			/*
			double kaon_weight = 1;
			double kaon_err = 0;
			if( applyKaon == 1 && this_bin_p > 1 ){
				kaon_weight = (double) kaonWeight[chargeIdx][this_bin_p]->GetBinContent(this_bin_xB, this_bin_Q2, this_bin_Z);
			}

			if( kaon_weight <= 0 || this_bin_p < 2 ){
				kaon_weight = 1;
			}
			*/
			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );

			double weight = 1;
			weight *= corrector.getCorrectionFactor(1, chargeIdx);
			if( weight > 5 ){continue;}
			weight *= corrector.getCorrectionFactor(2, chargeIdx);	

			hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ(), weight );
			hZ_uw[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ()  );
			//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			

		}
	}
	
	TTreeReader reader_k("ePi", kFile);

	TTreeReaderValue<electron> e_k(reader_k, "e");
	TTreeReaderArray<pion> k(reader_k, "pi");
	
	TTreeReaderArray<bool> isGoodKaon(reader_k, "isGoodPion");

	event_total = reader_k.GetEntries();
	//double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	while (reader_k.Next()) {
		int event_count = reader_k.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
		}

		for( int i = 0; i < (int) ( k.end() - k.begin() ); i++ ){
			
			int chargeIdx = (int)( k[i].getCharge() < 1 );
			double p_pi = k[i].get3Momentum().Mag();
			int this_bin_Q2 = (int)( ( (e_k->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e_k->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (k[i].getZ() - .3)/(1.-.3) )*bins_Z);

			bool matching = true;
			if( matchType == 2 ){ matching = !isGoodPion[i]; }

			if( matching ){ continue; }

			corrector.setKinematics( e_k->getXb(), e_k->getQ2(), k[i].getZ(), p_pi );

			double weight = 1;
			weight *= corrector.getCorrectionFactor(2, chargeIdx);
			if( weight <= 0 ){continue;}
			if( weight >= 1 ){ continue;}
			weight *= corrector.getCorrectionFactor(1, chargeIdx);	
			if( weight > 5 ){continue;}
			hZ_k_uw[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( k[i].getZ() );
			hZ_k[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( k[i].getZ(), weight );
			//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			

		}
	}

	outFile->cd();
	
	for( int i = 1; i <= bins_Q2; i++ ){
		for( int j = 1; j <= bins_xB; j++ ){

			hZ[i-1][j-1][0]->Divide( hZ[i-1][j-1][1] );
			hZ_uw[i-1][j-1][0]->Divide( hZ_uw[i-1][j-1][1] );
			hZ_k[i-1][j-1][0]->Divide( hZ_k[i-1][j-1][1] );
			hZ_k_uw[i-1][j-1][0]->Divide( hZ_k_uw[i-1][j-1][1] );

			hZ[i-1][j-1][0]->Write();
			//hZ[i-1][j-1][1]->Write();
			hZ_k[i-1][j-1][0]->Write();
			//hZ_k[i-1][j-1][1]->Write();
			
			hZ_uw[i-1][j-1][0]->Write();
			//hZ_uw[i-1][j-1][1]->Write();
			hZ_k_uw[i-1][j-1][0]->Write();
			//hZ_k_uw[i-1][j-1][1]->Write();
		}
	}		

	outFile->Close();

}
