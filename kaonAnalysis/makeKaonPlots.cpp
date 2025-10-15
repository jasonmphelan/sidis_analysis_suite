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
#include "TChain.h"
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
#include "analyzer.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"
#define HIST_PATH _HIST

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		cerr << "[Type (0 - pi2k, 1 - k2pi)] [Theta Cut (optional)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	int nBinsQ2 = bins_Q2;
	int nBinsXb = bins_xB/2;//(int) (10* (bins_xB/25.));
	int nBinsZ = bins_Z;

	//TString in_name = argv[1];
       	TString out_name = argv[1];
		
	int type = atoi(argv[2]);

       	double theta_cut = 0;
		if( argc > 3 ){ theta_cut = atoi(argv[3]); }

       	TFile * outFile = new TFile(out_name, "RECREATE");//(TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
	
	// Declare histograms

	cout<<"Creating Histograms\n";


	TH1F * h_Beta[2][nBinsQ2+1][nBinsXb+1][nBinsZ+1][bins_p+1];
	TH1F * h_Beta_rich[2][nBinsQ2+1][nBinsXb+1][nBinsZ+1][bins_p+1];
	TH2F * hBeta_p[2][nBinsQ2+1][nBinsXb+1][nBinsZ+1];
	TH2F * hBeta_rich_p[2][nBinsQ2+1][nBinsXb+1][nBinsZ+1];
	
	TString data_type[2] = {"pip", "pim"};

	for( int j = 0; j <= nBinsQ2; j++ ){
		for( int k = 0; k <= nBinsXb; k++ ){
			for( int l = 0; l <= nBinsZ; l++ ){
				for( int i = 0; i < 2; i++ ){//Bin by charge
					for( int m = 0; m <= bins_p; m++ ){//momentum bins
						h_Beta[i][j][k][l][m]             = new TH1F("hBeta_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_%i_%i;#beta;Counts [a.u.]", j, k), 135, -.1, 1.25);
						h_Beta_rich[i][j][k][l][m]             = new TH1F("hBeta_rich_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_rich_%i_%i;#beta;Counts [a.u.]", j, k), 135, -.1, 1.25);
					}
					hBeta_p[i][j][k][l]		= new TH2F("hBeta_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 500, .9 ,1.05 );
					hBeta_rich_p[i][j][k][l]		= new TH2F("hBeta_rich_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 500, 0.9 ,1.05 );
			
				}
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	//TFile * file = new TFile(in_name);
	TChain * file = new TChain("ePi");
	if(type == 0){
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root");
	}
	if(type == 1){
		cout<<"Adding kaon files"<<std::endl;
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.2/final_skim.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.4/final_skim.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.6/final_skim.root");
	}
		//TTreeReader reader("ePi", file);
	TTreeReader reader( file);

	TTreeReaderValue<electron> e(reader, "e");

	TTreeReaderArray<pion> pi(reader, "pi");
	
	TTreeReaderArray<bool> isGoodPion_vec(reader, "isGoodPion");
	//TTreeReaderArray<bool> isGoodPion_3d_vec(reader, "isGoodPion_3d");
        //TTreeReaderArray<bool> isGoodPion_no_acc_vec(reader, "isGoodPion_no_acc");

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);
	if( type == 0 ){
		anal.loadMatchingFunctions("matchCutPi2K.root");
	}
	if( type == 1 ){
		anal.loadMatchingFunctions("matchCutK2Pi.root");
	}
	int event_count = 0;
	while (reader.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
	       	event_count++;
                double Q2 = e->getQ2();
                double xB = e->getXb();


		int chargeIdx = 0;
	
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			//if(accMatchType < 2 && !isGoodPion_no_acc_vec[i]) {continue;}
			if( !isGoodPion_vec[i]) {continue;}
			//if(accMatchType == 3 && !isGoodPion_3d_vec[i]) {continue;}
			if(pi[i].getBeta_rich() < .0001){continue;}
			//if( pi[i].get3Momentum().Theta()*rad_to_deg > max_theta_cut ){ continue; }
			
			chargeIdx = (int)(pi[i].getCharge() < 1);
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta();
			double Z = pi[i].getZ();
			
			if( theta_cut > 0 && !anal.acceptance_match_2d(theta_pi*rad_to_deg, p_pi, 1) ){continue;}
			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nBinsXb) + 1;
			int this_bin_Z = (int)( ( (Z - Z_min)/(Z_max - Z_min) )*nBinsZ) + 1;
			int this_bin_p = -1;

			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
					this_bin_p = j+1;
				}
			}
		
			double m2_TOF = p_pi*p_pi*( 1 - pi[i].getBeta()*pi[i].getBeta() )/(pi[i].getBeta()*pi[i].getBeta());
			double m2_RICH = p_pi*p_pi*( 1 - pi[i].getBeta_rich()*pi[i].getBeta_rich() )/(pi[i].getBeta_rich()*pi[i].getBeta_rich());


			h_Beta[chargeIdx][0][0][0][0]->Fill(m2_TOF);
			h_Beta_rich[chargeIdx][0][0][0][0]->Fill(m2_RICH);
			hBeta_p[chargeIdx][0][0][0]->Fill(p_pi, pi[i].getBeta());
			hBeta_rich_p[chargeIdx][0][0][0]->Fill(p_pi, pi[i].getBeta_rich());
			
			h_Beta[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(m2_TOF);
			h_Beta_rich[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(m2_RICH);
			hBeta_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, pi[i].getBeta());
			hBeta_rich_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, pi[i].getBeta_rich());
			
		}
	}
	

	outFile->cd();
	
	cout<<"BEGIN WRITING HISTOGRAMS\n";

	for( int j = 0; j <= nBinsQ2; j++ ){
		for( int k = 0; k <= nBinsXb; k++ ){
			for( int l = 0; l <= nBinsZ; l++ ){
				for( int i = 0; i < 2; i++ ){//Bin by charge
					//hBeta_p[i][j][k][l]->Write();
					if(  hBeta_rich_p[i][j][k][l]->Integral() != 0){
						hBeta_rich_p[i][j][k][l]->Write();
						hBeta_p[i][j][k][l]->Write();
					}
					for( int m = 0; m <= 4; m++ ){//momentum bins
						if( j > 0 && k > 0 && l > 0 && m > 0 && h_Beta_rich[i][j][k][l][m]->Integral()!= 0){
							h_Beta_rich[i][j][k][l][m]->Write();
						}
			
					}
				}
			}	
		}
	
	}
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}
