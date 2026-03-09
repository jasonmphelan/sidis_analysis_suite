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
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "constants.h"
#include "cut_values.h"

#define HIST_PATH _HIST
#define CORR_PATH _DATA

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;


double fitF( double *p, double *par);

int main( int argc, char** argv){

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Correction File (no extension)] [Fit File (full path)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";


	//Read in files, constants
	TString in_name = argv[1];
	TString out_name = argv[2];
	TString fit_name = argv[3];
	double sigma = atof(argv[4]);

	TFile * outFile = new TFile(out_name, "RECREATE");
	TFile * outFile_fits = new TFile(fit_name, "RECREATE");
	TFile * inFile = new TFile(in_name);

	int nBinsQ2 = bins_Q2;
	int nBinsXb = bins_xB/2;
	int nBinsZ = bins_Z;
	
	//Load hists for fitting
	// Output only: summed over l
	TH2F* hBeta_sumZ[2][13][8] = {{{nullptr}}};

	for( int j = 0; j <= nBinsQ2; j++ ){
	for( int k = 0; k <= nBinsXb; k++ ){
		for( int l = 0; l <= nBinsZ; l++ ){

		// species 0: pip
		{
			TString name = TString("hBeta_rich_p_pip") + Form("_Q2_%i_xB_%i_Z_%i", j, k, l);
			TH2F* h_in = (TH2F*) inFile->Get(name);
			if( h_in ){
			// Make an owned copy so we're not tied to the TFile's directory
			TH2F* h_tmp = (TH2F*) h_in->Clone(Form("%s__tmp", name.Data()));
			h_tmp->SetDirectory(nullptr);

			if( !hBeta_sumZ[0][j][k] ){
				hBeta_sumZ[0][j][k] = (TH2F*) h_tmp->Clone(Form("hBeta_sumZ_pip_Q2_%i_xB_%i", j, k));
				hBeta_sumZ[0][j][k]->SetDirectory(nullptr);
				hBeta_sumZ[0][j][k]->Reset();
				// optional but often good if weights were used:
				// hBeta_sumZ[0][j][k]->Sumw2();
			}

			hBeta_sumZ[0][j][k]->Add(h_tmp);
			delete h_tmp;
			}
		}

		// species 1: pim
		{
			TString name = TString("hBeta_rich_p_pim") + Form("_Q2_%i_xB_%i_Z_%i", j, k, l);
			TH2F* h_in = (TH2F*) inFile->Get(name);
			if( h_in ){
			TH2F* h_tmp = (TH2F*) h_in->Clone(Form("%s__tmp", name.Data()));
			h_tmp->SetDirectory(nullptr);

			if( !hBeta_sumZ[1][j][k] ){
				hBeta_sumZ[1][j][k] = (TH2F*) h_tmp->Clone(Form("hBeta_sumZ_pim_Q2_%i_xB_%i", j, k));
				hBeta_sumZ[1][j][k]->SetDirectory(nullptr);
				hBeta_sumZ[1][j][k]->Reset();
				// optional:
				// hBeta_sumZ[1][j][k]->Sumw2();
			}

			hBeta_sumZ[1][j][k]->Add(h_tmp);
			delete h_tmp;
			}
		}

		}
	}
	}


	int bins_p = hBeta_sumZ[0][1][1]->GetXaxis()->GetNbins();

	//Declare Fit Pointers
	TF1 *fitPim[nBinsQ2][nBinsXb][bins_p];
	TF1 *fitPip[nBinsQ2][nBinsXb][bins_p];
	
	//Fit params
	double min[4] = { 0 };
	double max[4] = { .08, .08, .08,  .08 };
	double min_k[4] = { .985, .985, .985, .985 };

	//Store correction factor
	TH3F * kaonCorr_p;
	TH3F * kaonCorr_m;
	TH3F * kaonCorr_full;

	kaonCorr_p = new TH3F( Form("hKaonCorrP"),Form("hKaonCorrP"), nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, bins_p, 1.25, 5); 
	kaonCorr_m = new TH3F( Form("hKaonCorrM"),Form("hKaonCorrM") , nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, bins_p, 1.25, 5); 
	kaonCorr_full = new TH3F( Form("hKaonCorr"), "hKaonCorr", nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, bins_p, 1.25, 5); 


	cout<<"BEGIN FITTING\n";
	
	

	//Loop through bins for fitting
	outFile_fits->cd();
	for( int j = 0; j < nBinsQ2; j++ ){
		for( int k = 0; k < nBinsXb; k++ ){
			for( int m = 0; m < bins_p; m++ ){//momentum bins
				//If not filled, continue
				if(!hBeta_sumZ[0][j+1][k+1]){continue;}
				if(!hBeta_sumZ[1][j+1][k+1]){continue;}
				

				TH1F * hBeta_sumZ_px_pip =
						(TH1F*)hBeta_sumZ[0][j+1][k+1]->ProjectionY(
						Form("hBeta_sumZ_px_pip_Q2_%i_xB_%i_p_%i", j, k, m),
						m+1 ,m+1,   // y-bin range (inclusive), 1..NbinsY = all visible bins
						"e"                   // keep proper errors
						);


				TH1F * hBeta_sumZ_px_pim =
						(TH1F*)hBeta_sumZ[1][j+1][k+1]->ProjectionY(
						Form("hBeta_sumZ_px_pim_Q2_%i_xB_%i_p_%i", j, k, m),
						m+1, m+1,   // y-bin range (inclusive), 1..NbinsY = all visible bins
						"e"                   // keep proper errors
						);



				cout<<"Found histogram\n";

				if(hBeta_sumZ_px_pip->Integral()<=0){continue;}
				if(hBeta_sumZ_px_pim->Integral()<=0){continue;}

				hBeta_sumZ_px_pip->Print("all");


				hBeta_sumZ_px_pip->Write();
				hBeta_sumZ_px_pim->Write();
				//fit positive pions

				double fitProb_pip = 0;

				fitPip[j][k][m] = new TF1(Form("f_pip_%i_%i_%i", j+1, k+1, m+1), "gaus", min[m], max[m]);
				fitPip[j][k][m]->SetParameters( hBeta_sumZ_px_pip->GetMaximum(), .139*.139, .01 ); 
				hBeta_sumZ_px_pip->Fit( Form("f_pip_%i_%i_%i", j+1, k+1, m+1), "R" );
				fitProb_pip = fitPip[j][k][m]->GetProb();	

				//fit negative pions	
				double fitProb_pim = 0;
				fitPim[j][k][m] = new TF1(Form("f_pim_%i_%i_%i", j+1, k+1, m+1), "gaus", min[m], max[m]);
				fitPim[j][k][m]->SetParameters( hBeta_sumZ_px_pim->GetMaximum(), .139*.139, .01 ); 
				hBeta_sumZ_px_pim->Fit( Form("f_pim_%i_%i_%i", j+1, k+1, m+1), "R" );
				fitProb_pim = fitPim[j][k][m]->GetProb();	
				
				//Get Pip Counts for correction
				double pip_num = 0;
				double pip_fit_mean = fitPip[j][k][m]->GetParameter(1);
				double pip_fit_std = fitPip[j][k][m]->GetParameter(2);
				int pip_bin_min = 1;//hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pip_fit_mean - 2*pip_fit_std );
				int pip_bin_max = 1;
				
				if( pip_fit_mean < -0.05 || pip_fit_mean > 0.075 || fitProb_pip <= 0 || fitProb_pip >= 1 || hBeta_sumZ_px_pip->Integral() <= 150 ||
							pip_fit_mean + sigma*pip_fit_std > .2){	
					pip_bin_max= hBeta_sumZ_px_pip->FindBin( 0.1*sigma/2.5 );
				}
				else{
					pip_bin_max = hBeta_sumZ_px_pip->FindBin( pip_fit_mean + sigma*pip_fit_std );
					fitPip[j][k][m]->Write();	
				}
			

				for( int bin = 1; bin <= pip_bin_max; bin++ ){
					pip_num+=hBeta_sumZ_px_pip->GetBinContent(bin);
				}

				double kp_num = hBeta_sumZ_px_pip->Integral(0, hBeta_sumZ_px_pip->GetNbinsX()) - pip_num;
				double pip_err = (1./pow( kp_num+pip_num, 2))*sqrt( pip_num*kp_num*(pip_num+kp_num ) );
				
				//Get Pim Counts for correction
				
				double pim_num = 0;
				double pim_fit_mean = fitPim[j][k][m]->GetParameter(1);
				double pim_fit_std = fitPim[j][k][m]->GetParameter(2);
				int pim_bin_min = 1;//hBeta[1][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pim_fit_mean - 2*pim_fit_std );
				int pim_bin_max = 1;//hBeta[1][j+1][k+1][l+1][m+1]->GetBin( pim_fit_mean + 2*pim_fit_std );
				
				if( pim_fit_mean < 0 || pim_fit_mean > 0.05 || fitProb_pim <= 0 || fitProb_pim >= 1 || hBeta_sumZ_px_pim->Integral() <= 150||
				pim_fit_mean + sigma*pim_fit_std > .2){	
					pim_bin_max= hBeta_sumZ_px_pim->FindBin( 0.1 * sigma/2.5 );
				}
				else{
					pim_bin_max = hBeta_sumZ_px_pim->FindBin( pim_fit_mean + sigma*pim_fit_std );
					fitPim[j][k][m]->Write();	
				}

				for( int bin = 1; bin <= pim_bin_max; bin++ ){
					pim_num+=hBeta_sumZ_px_pim->GetBinContent(bin);
				}

				double km_num = hBeta_sumZ_px_pim->Integral(0, hBeta_sumZ_px_pim->GetNbinsX()) - pim_num;
				double pim_err = (1./pow( km_num+pim_num, 2))*sqrt( pim_num*km_num*(pim_num+km_num ) );
			

				if(pip_num/(pip_num+kp_num) >= 0 && pip_num/(pip_num+kp_num) < 1.  ){
					kaonCorr_p->SetBinContent(k+1, j+1, m+1, pip_num/(pip_num+kp_num));
					kaonCorr_p->SetBinError(k+1, j+1, m+1, pip_err);
				}
				else{
					kaonCorr_p->SetBinContent(k+1, j+1, m+1, 0);
					kaonCorr_p->SetBinError(k+1, j+1, m+1, 0);
				}
				if(pim_num/(pim_num+km_num) >= 0  && pim_num/(pim_num+km_num) < 1 ){
					kaonCorr_m->SetBinContent(k+1, j+1, m+1, pim_num/(pim_num+km_num));
					kaonCorr_m->SetBinError(k+1, j+1, m+1, pim_err);
				}
				else{
					kaonCorr_m->SetBinContent(k+1, j+1, m+1, -1);
					kaonCorr_m->SetBinError(k+1, j+1, m+1, 0);
				}
				
				//fitPim[j][k][l][m]->Write();	
				//fitPip[j][k][l][m]->Write();	
		
			}
			
		}
	}	
	outFile_fits->Close();
	outFile->cd();

	kaonCorr_full = (TH3F*)kaonCorr_p->Clone(Form("hKaonCorr")); 
	kaonCorr_full->Divide(kaonCorr_m);
	kaonCorr_p->Write();
	kaonCorr_m->Write();
	kaonCorr_full->Write();
	kaonCorr_full->Print();
	

	outFile->Close();
}

