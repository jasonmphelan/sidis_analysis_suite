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
		cerr << "./code [Output File (no extension)]\n";
		cerr << "[Type (0 - pi2k, 1 - k2pi)] [Theta Cut (optional)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	int nBinsQ2 = bins_Q2;
	int nBinsXb =(int) (10* (bins_xB/25.));
	int nBinsZ = 2*bins_Z;

	//TString in_name = argv[1];
       	TString out_name = argv[1];
		
	int type = atoi(argv[2]);

       	double max_theta_cut = 999;
	if( argc > 3 ){ max_theta_cut = atoi(argv[4]); }

       	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
	
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
						h_Beta[i][j][k][l][m]             = new TH1F("hBeta_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_%i_%i;#beta;Counts [a.u.]", j, k), 100, .97, 1.01);
						h_Beta_rich[i][j][k][l][m]             = new TH1F("hBeta_rich_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m), Form("Beta_rich_%i_%i;#beta;Counts [a.u.]", j, k), 100, .97, 1.01);
					}
					hBeta_p[i][j][k][l]		= new TH2F("hBeta_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 100, .95 ,1 );
					hBeta_rich_p[i][j][k][l]		= new TH2F("hBeta_rich_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 100, .95 ,1 );
			
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


	int event_count = 0;
	while (reader.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
	       	event_count++;
                double Q2 = e->getQ2();
                double xB = e->getXb();


		int chargeIdx = 0;
	
		//if( e->getW() < 19 || e->getV() < 19 ){continue;}


		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			//if(accMatchType < 2 && !isGoodPion_no_acc_vec[i]) {continue;}
			if( !isGoodPion_vec[i]) {continue;}
			//if(accMatchType == 3 && !isGoodPion_3d_vec[i]) {continue;}
			if(pi[i].getBeta_rich() < .0001){continue;}
			if( pi[i].get3Momentum().Theta()*rad_to_deg > max_theta_cut ){ continue; }
			chargeIdx = (int)(pi[i].getCharge() < 1);
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta();
			double Z = pi[i].getZ();

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nBinsXb) + 1;
			int this_bin_Z = (int)( ( (Z - Z_min)/(Z_max - Z_min) )*nBinsZ) + 1;
			int this_bin_p = -1;

			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
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
	
	//file->Close();

	outFile->cd();
	
	cout<<"BEGIN WRITING HISTOGRAMS\n";

	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int j = 0; j <= nBinsQ2; j++ ){
		for( int k = 0; k <= nBinsXb; k++ ){
			for( int l = 0; l <= nBinsZ; l++ ){
				for( int i = 0; i < 2; i++ ){//Bin by charge
					//hBeta_p[i][j][k][l]->Write();
					if( j > 0 && k > 0 && l > 0 && hBeta_rich_p[i][j][k][l]->Integral() != 0){

						hBeta_rich_p[i][j][k][l]->SetTitleSize(10);
						hBeta_rich_p[i][j][k][l]->SetTitle( Form( "%.1f<Q^{2}<%.1f & %.2f<x_{B}<%.2f & %.2f<Z<%.2f" , 
											2 + .5*(j-1), 2 + .5*j,
										     	.1 + .05*(k-1), .1 + .05*k,
											.3 + .05*(l - 1), .3 + .05*l));
						hBeta_rich_p[i][j][k][l]->Write();
						hBeta_rich_p[i][j][k][l]->Draw();
						//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");
						//canvas.Clear();
					}
					for( int m = 0; m <= 4; m++ ){//momentum bins
						if( j > 0 && k > 0 && l > 0 && m > 0 && h_Beta_rich[i][j][k][l][m]->Integral()!= 0){

							h_Beta_rich[i][j][k][l][m]->SetTitleSize(10);
							h_Beta_rich[i][j][k][l][m]->SetTitle( Form( "%.1f<Q^{2}<%.1f & %.2f<x_{B}<%.2f & %.2f<Z<%.2f & %.2f<p<%.2f" , 
												2 + .5*(j-1), 2 + .5*j,
											     	.1 + .05*(k-1), .1 + .05*k,
												.3 + .05*(l - 1), .3 + .05*l,
												p_bin_edges[m-1], p_bin_edges[m]));
						h_Beta_rich[i][j][k][l][m]->Write();
						h_Beta_rich[i][j][k][l][m]->Draw();
						h_Beta[i][j][k][l][m]->SetLineColor(kRed);
					//	h_Beta[i][j][k][l][m]->Scale( (double)h_Beta_rich[i][j][k][l][m]->Integral()/(double)h_Beta[i][j][k][l][m]->Integral());
						h_Beta[i][j][k][l][m]->Draw("same");
						canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");
						canvas.Clear();
						}
			
					}
				}
			}	
		}
	
	}
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}
