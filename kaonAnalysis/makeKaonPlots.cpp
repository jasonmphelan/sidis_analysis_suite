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

#define CORR_PATH _DATA

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

	int accMatchType = 2;
	if( argc > 3 ){ accMatchType = atoi(argv[3]); }
	//TString corrFileName = "corrections.root";
	//if( argc > 4 ){ corrFileName = argv[4]; }
       
       	TFile * outFile = new TFile(out_name, "RECREATE");
	
	// Declare histograms

	cout<<"Creating Histograms\n";


	TH1F * h_Beta[2][bins_Q2+1][bins_xB+1][bins_Z+1][bins_p+1];
	TH1F * h_Beta_rich[2][bins_Q2+1][bins_xB+1][bins_Z+1][bins_p+1];
	TH2F * hBeta_p[2][bins_Q2+1][bins_xB+1][bins_Z+1];
	TH2F * hBeta_rich_p[2][bins_Q2+1][bins_xB+1][bins_Z+1];
	
	TString data_type[2] = {"pip", "pim"};

	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){
			for( int l = 0; l <= bins_Z; k++ ){
				for( int i = 0; i < 2; i++ ){//Bin by charge
					for( int m = 0; m <= 4; k++ ){//momentum bins
						h_Beta[i][j][k][l][m]             = new TH1F("hBeta_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_%i_%i;#beta;Counts [a.u.]", j, k), 50, .95, 1);
						h_Beta_rich[i][j][k][l][m]             = new TH1F("hBeta_rich_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_rich_%i_%i;#beta;Counts [a.u.]", j, k), 50, .95, 1);
					}
					hBeta_p[i][j][k][l]		= new TH2F("hBeta_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 100, .95 ,1 );
					hBeta_rich_p[i][j][k][l]		= new TH2F("hBeta_rich_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 100, .95 ,1 );
			
				}
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	TFile * file = new TFile(in_name);
	TTreeReader reader("ePi", file);

	TTreeReaderValue<electron> e(reader, "e");

	TTreeReaderArray<pion> pi(reader, "pi");
	
	TTreeReaderArray<bool> isGoodPion_vec(reader, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion_3d_vec(reader, "isGoodPion_3d");
        TTreeReaderArray<bool> isGoodPion_no_acc_vec(reader, "isGoodPion_no_acc");


	int event_count = 0;
	while (reader.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
                event_count++;

                double Q2 = e->getQ2();
                double xB = e->getXb();


		int chargeIdx = 0;
		
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			if(accMatchType < 2 && !isGoodPion_no_acc_vec[i]) {continue;}
			if(accMatchType == 2 && !isGoodPion_vec[i]) {continue;}
			if(accMatchType == 3 && !isGoodPion_3d_vec[i]) {continue;}
			
			chargeIdx = (int)(pi[i].getCharge() < 1);
			//double M_x = M_x_vec[i]; 
			//double pT_pi = pi_q_vec[i].Vect().Pt();
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta();
			double Z = pi[i].getZ();

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
			int this_bin_Z = (int)( ( (Z - Z_min)/(Z_max - Z_min) )*bins_Z) + 1;
			int this_bin_p = -1;

			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi > p_bin_edges[j+1] ){
					this_bin_p = j+1;
				}
			}
		
			h_Beta[chargeIdx][0][0][0][0]->Fill(pi[i].getBeta());
			h_Beta_rich[chargeIdx][0][0][0][0]->Fill(pi[i].getBeta_rich());
			hBeta_p[chargeIdx][0][0][0]->Fill(p_pi, pi[i].getBeta());
			hBeta_rich_p[chargeIdx][0][0][0]->Fill(p_pi, pi[i].getBeta_rich());
			
			h_Beta[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(pi[i].getBeta());
			h_Beta_rich[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(pi[i].getBeta_rich());
			hBeta_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, pi[i].getBeta());
			hBeta_rich_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, pi[i].getBeta_rich());
			
		}
	}
	
	file->Close();

	outFile->cd();
	
	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){
			for( int l = 0; l <= bins_Z; k++ ){
				for( int i = 0; i < 2; i++ ){//Bin by charge
					hBeta_p[i][j][k][l]->Write();
					hBeta_rich_p[i][j][k][l]->Write();
					for( int m = 0; m <= 4; k++ ){//momentum bins
						h_Beta[i][j][k][l][m]->Write();
						h_Beta_rich[i][j][k][l][m]->Write();
			
					}
				}
			}	
		}
	
	}
	outFile->Close();
}
