
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
		cerr << "[Acc/Bin Correction File]\n";
		cerr << "[Kaon Correction File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
	int matchType = atoi(argv[3]);
	TString corrFileName = argv[4];
	TString kaonFileName = argv[5];

       	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TFile * inFile = new TFile( in_name );
	TFile * weightFile = new TFile( "../corrections/corrections/" + corrFileName );	

	TH3F * accWeight_pip = (TH3F *)weightFile->Get("hAccCorrectionP");
	TH3F * accWeight_pim = (TH3F *)weightFile->Get("hAccCorrectionM");
	TH3F * binWeight_pip = (TH3F *)weightFile->Get("hBinMigrationP");
	TH3F * binWeight_pim = (TH3F *)weightFile->Get("hBinMigrationM");

	TFile * kaonFile = new TFile((TString) CORR_PATH + "/correctionFiles/corrections_kaons.root"+kaonFileName);
	TH3F * kaonWeight_pip[bins_p];
	TH3F * kaonWeight_pim[bins_p];

	for(int i = 0; i < bins_p; i++){
		kaonWeight_pip[i] = (TH3F *)weightFile->Get(Form( "hKaonCorrectionsP_%i", i ));
		kaonWeight_pim[i] = (TH3F *)weightFile->Get(Form( "hKaonCorrectionsM_%i", i ));
	}

	TH1F * hZ[bins_Q2+1][bins_xB+1][2];

	TString charge_str[2] = {"", "Pim"};

	for( int i = 0; i <= bins_Q2; i++ ){
		for( int j = 0; j <= bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i, j) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i, j) , bins_Z, .3, 1);
						
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
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z) + 1;
			int this_bin_p = -1;

			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
					this_bin_p = j+1;
				}
			}
			
			if( this_bin_p < 1 ){ continue; } //just in case

			double acc_weight = 1.;
			double bin_weight = 1.;
			double kaon_weight = 1.;
			
			double acc_err = 0.;
			double bin_err = 0.;
			double kaon_err = 0.;
			
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[i]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[i]; }
			//else{ matching = false; }

			if( matching ){ continue; }

			hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ() );
			events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			
			//hZ[0][0][chargeIdx]->Fill(Z, acc_weight*bin_weight);

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
			
			hZ[i][j][0]->Sumw2();
			hZ[i][j][1]->Sumw2();

				for( int k = 1; k <= bins_Z; k++ ){
					double n_pip_corr = 0;
					double pip_err = 0;
					double acc_weight_pip = (double) accWeight_pip->GetBinContent( j, i, k );
					double bin_weight_pip = (double) binWeight_pip->GetBinContent( j, i, k );
					double kaon_weight_pip;

					double acc_err_pip = (double) accWeight_pip->GetBinError( j, i, k );
					double bin_err_pip = (double) binWeight_pip->GetBinError( j, i, k );
					double kaon_err_pip;
				
					double err_term_1 = 0;
					double err_term_2 = 0;
					double err_term_3 = 0;

					for( int p = 0; p < bins_p; p++ ){
						kaon_weight_pip = (double) kaonWeight_pip[p]->GetBinContent(i, j, k);
						kaon_err_pip = (double) kaonWeight_pip[p]->GetBinError(i, j, k);
					
						if( kaon_weight_pip == 0 ){
							kaon_weight_pip = 1;
							kaon_err_pip = 0;
						}

						err_term_1 += kaon_weight_pip*events_in_bin[0][i][j][k][p];
						err_term_2 += pow( kaon_err_pip*events_in_bin[0][i][j][k][p], 2);
					     	err_term_3 += pow( kaon_weight_pip*sqrt(events_in_bin[0][i][j][k][p]), 2);	
					}

					n_pip_corr = acc_weight_pip*bin_weight_pip*err_term_1;
					pip_err = sqrt( pow( acc_weight_pip*bin_err_pip*err_term_1, 2) 
							+ pow(bin_weight_pip*acc_err_pip*err_term_1, 2)
							+ pow(bin_weight_pip*acc_weight_pip, 2)*err_term_2
							+ pow(bin_weight_pip*acc_weight_pip, 2)*err_term_3 );
					
					double n_pim_corr = 0;
					double pim_err = 0;
					double acc_weight_pim = (double) accWeight_pim->GetBinContent( j, i, k );
					double bin_weight_pim = (double) binWeight_pim->GetBinContent( j, i, k );
					double kaon_weight_pim;

					double acc_err_pim = (double) accWeight_pim->GetBinError( j, i, k );
					double bin_err_pim = (double) binWeight_pim->GetBinError( j, i, k );
					double kaon_err_pim;
				
					err_term_1 = 0;
					err_term_2 = 0;
					err_term_3 = 0;

					for( int p = 0; p < bins_p; p++ ){
						kaon_weight_pim = (double) kaonWeight_pim[p]->GetBinContent(i, j, k);
						kaon_err_pim = (double) kaonWeight_pim[p]->GetBinError(i, j, k);
					
						if( kaon_weight_pim == 0 ){
							kaon_weight_pim = 1;
							kaon_err_pim = 0;
						}

						err_term_1 += kaon_weight_pim*events_in_bin[1][i][j][k][p];
						err_term_2 += pow( kaon_err_pim*events_in_bin[1][i][j][k][p], 2);
					     	err_term_3 += pow( kaon_weight_pim*sqrt(events_in_bin[1][i][j][k][p]), 2);	
					}

					n_pim_corr = acc_weight_pim*bin_weight_pim*err_term_1;
					pim_err = sqrt( pow( acc_weight_pim*bin_err_pim*err_term_1, 2) 
							+ pow(bin_weight_pim*acc_err_pim*err_term_1, 2)
							+ pow(bin_weight_pim*acc_weight_pim, 2)*err_term_2
							+ pow(bin_weight_pim*acc_weight_pim, 2)*err_term_3 );
							


					/*
					double n_pip = hZ[i][j][0]->GetBinContent( k );
					double n_pim = hZ[i][j][1]->GetBinContent( k );

						
					double n_pip_err = sqrt(n_pip);			
					double n_pim_err = sqrt(n_pim);			

					double acc_weight_pip = (double) accWeight_pip->GetBinContent( j, i, k );
					double bin_weight_pip = (double) binWeight_pip->GetBinContent( j, i, k );

					
					double acc_err_pip = (double) accWeight_pip->GetBinError( j, i, k );
					double bin_err_pip = (double) binWeight_pip->GetBinError( j, i, k );
					
					double acc_weight_pim = (double) accWeight_pim->GetBinContent( j, i, k );
					double acc_err_pim = (double) accWeight_pim->GetBinError( j, i, k );

					double bin_err_pim = (double) binWeight_pim->GetBinError( j, i, k );
					double bin_weight_pim = (double) binWeight_pim->GetBinContent( j, i, k );

					
					if( !isfinite(acc_weight_pip) || acc_err_pip/acc_weight_pip > .2){acc_weight_pip =0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
					if( !isfinite(bin_weight_pip) || bin_err_pip/bin_weight_pip > .2){bin_weight_pip= 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
					
					if( !isfinite(acc_weight_pim) || acc_err_pim/acc_weight_pim > .2){acc_weight_pim =0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
					if( !isfinite(bin_weight_pim) || bin_err_pim/bin_weight_pim > .2){bin_weight_pim= 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
				
					double n_pip_corr = n_pip*acc_weight_pip*bin_weight_pip;
					double pip_err = n_pip_corr*sqrt( pow( 1./sqrt(n_pip) , 2 ) + pow( acc_err_pip/acc_weight_pip, 2 ) + pow(bin_err_pip/bin_weight_pip, 2 ) );	
		
					double n_pim_corr = n_pim*acc_weight_pim*bin_weight_pim;
					double pim_err = n_pim_corr*sqrt( pow( 1./sqrt(n_pim) , 2 ) + pow( acc_err_pim/acc_weight_pim, 2 ) + pow(bin_err_pim/bin_weight_pim, 2 ) );	
					*/
					
					hZ[i][j][0]->SetBinContent(k, n_pip_corr);
					hZ[i][j][1]->SetBinContent(k, n_pim_corr);
					
					hZ[i][j][0]->SetBinError(k, pip_err);
					hZ[i][j][1]->SetBinError(k, pim_err);
				
						

				}


			hZ[i][j][0]->Divide(hZ[i][j][1]);

			TH1F * helper_3 = (TH1F *)hZ[i][j][0]->Clone();
			hZ[i][j][0]->Scale(-1.);
			hZ[i][j][0]->Add(helper_2);

			helper_3->Scale(4.);
			helper_3->Add( helper_1 );

			hZ[i][j][0]->Divide(helper_3);
			
			//for( int k = 1; k <= bins_Z; k++ ){

			//	if( isnan(hZ[i][j][0]->GetBinContent(k))

			hZ[i][j][0]->Write();
			hZ[i][j][1]->Write();
		}
	}		
	
	outFile->Close();
	
}
