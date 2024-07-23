#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"


#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<clashit>+;
#pragma link C++ class vector<clashit>+;
#endif

using std::cerr;
using std::isfinite;
using std::isinf;
using std::cout;
using std::ofstream;

const int nXbBins = 10;
const int nQ2Bins = 12;
const int nZBins = 14;

const double xB_min = .1;
const double xB_max = .6;
const double xB_width = (xB_max - xB_min)/( (double) nXbBins);

const double Q2_min = 2;
const double Q2_max = 8;
const double Q2_width = (Q2_max - Q2_min)/( (double) nQ2Bins);


void makeRatioBinned(TString corrFileName, TString inName, TString outName, int weights, int matchType){

	TFile * outFile = new TFile( outName, "RECREATE");
	TFile * inFile = new TFile( inName );
	TFile * weightFile = new TFile( "../corrections/corrections/" + corrFileName );	

	TH3F * accWeight_pip = (TH3F *)weightFile->Get("hAccCorrectionP");
	TH3F * accWeight_pim = (TH3F *)weightFile->Get("hAccCorrectionM");
	TH3F * binWeight_pip = (TH3F *)weightFile->Get("hBinMigrationP");
	TH3F * binWeight_pim = (TH3F *)weightFile->Get("hBinMigrationM");

	TH1F * hZ[nQ2Bins+1][nXbBins+1][2];

	TString charge_str[2] = {"", "Pim"};

	for( int i = 0; i <= nQ2Bins; i++ ){
		for( int j = 0; j <= nXbBins; j++ ){
			for( int k = 0; k < 2; k++ ){
		
				hZ[i][j][k] = new TH1F("hRatio" + charge_str[k] + Form("_%i_%i",  i, j) ,"hRatio_" + charge_str[k] + Form("_%i_%i",  i, j) , 14, .3, 1);
						
			}
		}
	}

	TTreeReader reader_rec("ePi", inFile);

	TTreeReaderValue<double> Q2_ptr(reader_rec, "Q2");
	TTreeReaderValue<double> xB_ptr(reader_rec, "xB");

	TTreeReaderArray<double> Z_vec(reader_rec, "Z");
	TTreeReaderArray<int> charge(reader_rec, "charge");
	
	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");



	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<endl;
		}


		//if( nLeadIdx < 0 ) continue;

		double Q2 = *Q2_ptr;
		double xB = *xB_ptr;
	

		for( int i = 0; i < (int) ( Z_vec.end() - Z_vec.begin() ); i++ ){
			double Z = Z_vec[i];

			int chargeIdx = (int)( charge[i] < 1 );

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins) + 1;
			int this_bin_Z = (int)( ( (Z - .3)/(1.-.3) )*nZBins) + 1;
		
			double acc_weight = 1.;
			double bin_weight = 1.;
			
			double acc_err = 0.;
			double bin_err = 0.;
			
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[i]; }
			else if( matchType == 3 ){ matching = !isGoodPion3d[i]; }
			else{ matching = false; }

			if( matching ){ continue; }

			hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill(Z);
			//hZ[0][0][chargeIdx]->Fill(Z, acc_weight*bin_weight);

		}
	}

	outFile->cd();
	TH1F * helper_1 = new TH1F("helper1", "helper1", 14, .3, 1);
	TH1F * helper_2 = new TH1F("helper2", "helper2", 14, .3, 1);
	
	for( int i = 1; i <= 14; i++ ){
		helper_1->SetBinContent(i, -1.);
		helper_2->SetBinContent(i, 4.);
		
		helper_1->SetBinError(i, 0.);
		helper_2->SetBinError(i, 0.);
	}

	for( int i = 1; i <= nQ2Bins; i++ ){
		for( int j = 1; j <= nXbBins; j++ ){
			
			hZ[i][j][0]->Sumw2();
			hZ[i][j][1]->Sumw2();

			if( weights == 1 ){
				for( int k = 1; k <= nZBins; k++ ){
				
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
				
					
					hZ[i][j][0]->SetBinContent(k, n_pip_corr);
					hZ[i][j][1]->SetBinContent(k, n_pim_corr);
					
					hZ[i][j][0]->SetBinError(k, pip_err);
					hZ[i][j][1]->SetBinError(k, pim_err);
				
						

				}

			}

			hZ[i][j][0]->Divide(hZ[i][j][1]);

			TH1F * helper_3 = (TH1F *)hZ[i][j][0]->Clone();
			hZ[i][j][0]->Scale(-1.);
			hZ[i][j][0]->Add(helper_2);

			helper_3->Scale(4.);
			helper_3->Add( helper_1 );

			hZ[i][j][0]->Divide(helper_3);
			
			//for( int k = 1; k <= nZBins; k++ ){

			//	if( isnan(hZ[i][j][0]->GetBinContent(k))

			hZ[i][j][0]->Write();
			hZ[i][j][1]->Write();
		}
	}		
	
	outFile->Close();
	
}
