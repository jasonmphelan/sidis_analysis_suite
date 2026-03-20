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
#include "DCfid_SIDIS.h"
#define HIST_PATH _HIST

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut );


int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		cerr << "[Type (0 - pi2k, 1 - k2pi)] [Theta Cut (optional)] [Mx_2pi^2 cut (optional, default 1.25)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	int nBinsQ2 = bins_Q2;
	int nBinsXb = bins_xB/2;//(int) (10* (bins_xB/25.));
	int nBinsZ = bins_Z;

    TString out_name = argv[1];
		
	int type = atoi(argv[2]);

    double theta_cut = 0;
	if( argc > 3 ){ theta_cut = atoi(argv[3]); }
	double Mx_2pi_cut = 1.25;
	if( argc > 4 ){ Mx_2pi_cut = atof(argv[4]); }

	// Optional: specify Q2 and xB bin ranges to sum over (1-indexed, inclusive)
	int Q2_bin_lo = 1, Q2_bin_hi = nBinsQ2;
	int xB_bin_lo = 1, xB_bin_hi = nBinsXb;
	if( argc > 5 ){ Q2_bin_lo = std::max(1, std::min(atoi(argv[5]), nBinsQ2)); }
	if( argc > 6 ){ Q2_bin_hi = std::max(1, std::min(atoi(argv[6]), nBinsQ2)); }
	if( argc > 7 ){ xB_bin_lo = std::max(1, std::min(atoi(argv[7]), nBinsXb)); }
	if( argc > 8 ){ xB_bin_hi = std::max(1, std::min(atoi(argv[8]), nBinsXb)); }
	cerr << Form("Summing Q2 bins [%d-%d], xB bins [%d-%d]\n", Q2_bin_lo, Q2_bin_hi, xB_bin_lo, xB_bin_hi);

	int target = 0;  // 0 = RGB/deuterium, 1 = RGA/proton
	if( argc > 9 ){ target = atoi(argv[9]); }
	cerr << "Target: " << target << " (" << (target == 1 ? "RGA/proton" : "RGB/deuterium") << ")\n";


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
					hBeta_p[i][j][k][l]		= new TH2F("hBeta_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 100, 1.25, 5, 135, -.1 ,1.25 );
					hBeta_rich_p[i][j][k][l]		= new TH2F("hBeta_rich_p_"+data_type[i]+Form("_Q2_%i_xB_%i_Z_%i", j, k, l), "", 20, 3, 5, 135, -.1 ,1.25 );
			
				}
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	//TFile * file = new TFile(in_name);
	TChain * file = new TChain("ePi");
	if(type == 0){
		if( target == 0 ){
			file->Add("../trees/final_skims/10.2/tight_skim.root");
			file->Add("../trees/final_skims/10.4/tight_skim.root");
			file->Add("../trees/final_skims/10.6/tight_skim.root");
		}
		if( target == 1 ){
			file->Add("../trees/final_skims/10.6/final_skim_rga.root");
		}
	}
	if(type == 1){
		cout<<"Adding kaon files"<<std::endl;
		if( target == 0 ){
			file->Add("../trees/final_skims/kaons_10.2/tight_skim.root");
			file->Add("../trees/final_skims/kaons_10.4/tight_skim.root");
			file->Add("../trees/final_skims/kaons_10.6/tight_skim.root");
		}
		if( target == 1 ){
			file->Add("../trees/final_skims/kaons_10.6/final_skim_rga.root");
		}
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

	anal.setAnalyzerLevel(0);//runType);

	if( type == 0 ){
		anal.loadMatchingFunctions("matchCutPi2K_map.root");
	}
	if( type == 1 ){
		anal.loadMatchingFunctions("matchCutK2Pi_map.root");
	}
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));
	int event_count = 0;
	while (reader.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
	       	event_count++;
                double Q2 = e->getQ2();
                double xB = e->getXb();


		int chargeIdx = 0;
		if( type == 0 &&  (int) ( pi.end() - pi.begin() ) == 2
				 && Mx_2pi_cut > 0 && 
				 checkDiffractive(pi[0], pi[1], *e, Mx_2pi_cut)) continue;

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			if( !isGoodPion_vec[i]) {continue;}
			if(pi[i].getBeta_rich() < .0001){continue;}
			
			chargeIdx = (int)(pi[i].getCharge() < 1);
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta();
			double phi_pi = pi[i].get3Momentum().Phi();
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
			hBeta_rich_p[chargeIdx][0][0][0]->Fill(p_pi, m2_RICH);
			
			h_Beta[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(m2_TOF);
			h_Beta_rich[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]->Fill(m2_RICH);
			hBeta_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, pi[i].getBeta());
			hBeta_rich_p[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z]->Fill(p_pi, m2_RICH);
			
		}
	}


	outFile->cd();
	
	// Build summed histograms over the specified Q2/xB bin ranges
	cout << "Creating summed histograms\n";
	TString sum_label = Form("_Q2_%d-%d_xB_%d-%d", Q2_bin_lo, Q2_bin_hi, xB_bin_lo, xB_bin_hi);
	TH1F * h_Beta_sum[2][nBinsZ+1][bins_p+1];
	TH1F * h_Beta_rich_sum[2][nBinsZ+1][bins_p+1];
	TH2F * hBeta_p_sum[2][nBinsZ+1];
	TH2F * hBeta_rich_p_sum[2][nBinsZ+1];
	for( int i = 0; i < 2; i++ ){
		for( int l = 0; l <= nBinsZ; l++ ){
			hBeta_p_sum[i][l]      = new TH2F("hBeta_p_sum_"+data_type[i]+Form("_Z_%i",l)+sum_label, "", 100, 1.25, 5, 135, -.1, 1.25);
			hBeta_rich_p_sum[i][l] = new TH2F("hBeta_rich_p_sum_"+data_type[i]+Form("_Z_%i",l)+sum_label, "", 20, 3, 5, 135, -.1, 1.25);
			for( int m = 0; m <= bins_p; m++ ){
				h_Beta_sum[i][l][m]      = new TH1F("hBeta_sum_"+data_type[i]+Form("_Z_%i_p_%i",l,m)+sum_label, "m^{2};Counts [a.u.]", 135, -.1, 1.25);
				h_Beta_rich_sum[i][l][m] = new TH1F("hBeta_rich_sum_"+data_type[i]+Form("_Z_%i_p_%i",l,m)+sum_label, "m^{2};Counts [a.u.]", 135, -.1, 1.25);
			}
		}
	}
	for( int j = Q2_bin_lo; j <= Q2_bin_hi; j++ ){
		for( int k = xB_bin_lo; k <= xB_bin_hi; k++ ){
			for( int i = 0; i < 2; i++ ){
				for( int l = 0; l <= nBinsZ; l++ ){
					hBeta_p_sum[i][l]->Add(hBeta_p[i][j][k][l]);
					hBeta_rich_p_sum[i][l]->Add(hBeta_rich_p[i][j][k][l]);
					for( int m = 0; m <= bins_p; m++ ){
						h_Beta_sum[i][l][m]->Add(h_Beta[i][j][k][l][m]);
						h_Beta_rich_sum[i][l][m]->Add(h_Beta_rich[i][j][k][l][m]);
					}
				}
			}
		}
	}

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
	// Write summed histograms
	for( int i = 0; i < 2; i++ ){
		for( int l = 0; l <= nBinsZ; l++ ){
			if( hBeta_rich_p_sum[i][l]->Integral() != 0 ){
				hBeta_p_sum[i][l]->Write();
				hBeta_rich_p_sum[i][l]->Write();
			}
			for( int m = 1; m <= bins_p; m++ ){
				if( l > 0 && h_Beta_rich_sum[i][l][m]->Integral() != 0 ){
					h_Beta_sum[i][l][m]->Write();
					h_Beta_rich_sum[i][l][m]->Write();
				}
			}
		}
	}
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}


bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut ){
	TLorentzVector missing = p_rest + e.getQ() - pi_1.get4Momentum() - pi_2.get4Momentum();
	return missing.Mag2() < Mx_cut;
}