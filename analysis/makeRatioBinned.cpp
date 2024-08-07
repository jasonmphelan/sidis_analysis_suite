
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
		cerr << "./code [Input File] [Output File]\n";
		cerr << "[Acceptance Matching Type (2,3 etc)] \n";
		cerr << "[Acc/Bin Correction?]\n";
		cerr << "[Acc/Bin Correction File]\n";
		cerr << "[Kaon Correction?]\n";
		cerr << "[Kaon Correction File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
	int matchType = atoi(argv[3]);
	int applyCorr = atoi(argv[4]);
	int applyKaon;
	TString corrFileName;
	TString kaonFileName;
	if(applyCorr == 1){
		corrFileName = argv[5];
		applyKaon = atoi( argv[6] );
	}
	else{ applyKaon = atoi( argv[5] ); }
	if(applyKaon == 1 && applyCorr == 1){
		kaonFileName = argv[7];
	}
	else if( applyKaon == 1 && applyCorr == 0 ){
		kaonFileName = argv[6];
	}

       	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TFile * inFile = new TFile( in_name );
	TFile * weightFile;
	TFile * kaonFile;

	TH3F * accWeight[2];
	TH3F * binWeight[2];
	TH3F * kaonWeight[2][bins_p];

	if( applyCorr == 1 ){
		weightFile = new TFile( (TString) CORR_PATH + "/correctionFiles/" + corrFileName );	
		accWeight[0] = (TH3F *)weightFile->Get("hAccCorrectionP");
		accWeight[1] = (TH3F *)weightFile->Get("hAccCorrectionM");
		binWeight[0] = (TH3F *)weightFile->Get("hBinMigrationP");
		binWeight[1] = (TH3F *)weightFile->Get("hBinMigrationM");
	}

	if( applyKaon == 1 ){
		kaonFile = new TFile((TString) CORR_PATH + "/correctionFiles/"+kaonFileName);
		for(int i = 0; i < bins_p; i++){
			kaonWeight[0][i] = (TH3F *)kaonFile->Get(Form( "hKaonCorrP_%i", i ));
			kaonWeight[1][i] = (TH3F *)kaonFile->Get(Form( "hKaonCorrM_%i", i ));
		}
	}

	TH1F * hZ[bins_Q2+1][bins_xB+1][2];

	TString charge_str[2] = {"", "Pim"};

	for( int i = 0; i < bins_Q2; i++ ){
		for( int j = 0; j < bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i+1, j+1) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , bins_Z, .3, 1);
						
			}
		}
	}

	TTreeReader reader_rec("ePi", inFile);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	
	//TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");
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
			int this_bin_p = -1;

			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
					this_bin_p = j;
				}
			}
			
			if( this_bin_p < 0 ){ continue; } //just in case

			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[i]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[i]; }
			//else{ matching = false; }

			if( matching ){ continue; }

			//hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ() );
			events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			

		}
	}

	outFile->cd();
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
			
			cout<<"------------------------------------------------\n";
			cout<<"Q2 Bin : "<<i<<std::endl;
			cout<<"xB Bin : "<<j<<std::endl;
			cout<<"------------------------------------------------\n\n";
			for( int k = 1; k <= bins_Z; k++ ){
				for( int charge = 0; charge <= 1; charge++ ){
					double n_pi_corr = 0;
					double pi_err = 0;
					double acc_weight = 1;
					double bin_weight = 1;
					double acc_err = 0;
					double bin_err = 0;
				
					if( applyCorr == 1 ){
						acc_weight = (double) accWeight[charge]->GetBinContent( j, i, k );
						bin_weight = (double) binWeight[charge]->GetBinContent( j, i, k );
						bin_err = (double) binWeight[charge]->GetBinError( j, i, k );
						acc_err = (double) accWeight[charge]->GetBinError( j, i, k );
						
						if( !isfinite(acc_weight) || acc_err/acc_weight > .2){acc_weight = 0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
						if( !isfinite(bin_weight) || bin_err/bin_weight > .2){bin_weight = 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
					}

					double kaon_weight = 1;
					double kaon_err = 0;
				
					double err_term_1 = 0;
					double err_term_2 = 0;
					double err_term_3 = 0;

					for( int p = 0; p < bins_p; p++ ){
						if( applyKaon == 1 ){
							kaon_weight = (double) kaonWeight[charge][p]->GetBinContent(j, i, k);
							kaon_err = (double) kaonWeight[charge][p]->GetBinError(j, i, k);
						}

						if( kaon_weight == 0 ){
							kaon_weight = 1;
							kaon_err = 0;
						}

						err_term_1 += kaon_weight*events_in_bin[charge][i-1][j-1][k-1][p];
						err_term_2 += pow( kaon_err*events_in_bin[charge][i-1][j-1][k-1][p], 2);
						err_term_3 += pow( kaon_weight*sqrt(events_in_bin[charge][i-1][j-1][k-1][p]), 2);	
					}

					n_pi_corr = acc_weight*bin_weight*err_term_1;
					pi_err = sqrt( pow( acc_weight*bin_err*err_term_1, 2) 
							+ pow(bin_weight*acc_err*err_term_1, 2)
							+ pow(bin_weight*acc_weight, 2)*err_term_2
							+ pow(bin_weight*acc_weight, 2)*err_term_3 );

					hZ[i-1][j-1][charge]->SetBinContent(k, n_pi_corr);
					hZ[i-1][j-1][charge]->SetBinError(k, pi_err);
							

				}
			
			}

			hZ[i-1][j-1][0]->Divide(hZ[i-1][j-1][1]);

			TH1F * helper_3 = (TH1F *)hZ[i-1][j-1][0]->Clone();
			hZ[i-1][j-1][0]->Scale(-1.);
			hZ[i-1][j-1][0]->Add(helper_2);

			helper_3->Scale(4.);
			helper_3->Add( helper_1 );

			hZ[i-1][j-1][0]->Divide(helper_3);
		
			hZ[i-1][j-1][0]->Print("ALL");
			hZ[i-1][j-1][0]->Write();
			hZ[i-1][j-1][1]->Write();
		}
	}		

	outFile->Close();

}
