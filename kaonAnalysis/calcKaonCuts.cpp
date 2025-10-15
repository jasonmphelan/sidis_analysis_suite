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
	TH1F * hBeta[2][nBinsQ2+1][nBinsXb+1][nBinsZ+1][bins_p+1];

	for( int j = 0; j <= nBinsQ2; j++ ){
		for( int k = 0; k <= nBinsXb; k++ ){
			for( int l = 0; l <= nBinsZ; l++ ){
				for( int m = 0; m <= bins_p; m++ ){//momentum bins
					hBeta[0][j][k][l][m] = (TH1F*)inFile->Get((TString)"hBeta_rich_pip"+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m)); 
					hBeta[1][j][k][l][m] = (TH1F*)inFile->Get((TString)"hBeta_rich_pim"+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m)); 
					
				}
			}
		}
	}

	//Declare Fit Pointers
	TF1 *fitPim[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	TF1 *fitPip[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	
	//Fit params
	double min[4] = { 0 };
	double max[4] = { .08, .08, .08,  .08 };
	double min_k[4] = { .985, .985, .985, .985 };

	//Store correction factor
	TH3F * kaonCorr_p[bins_p];
	TH3F * kaonCorr_m[bins_p];
	TH3F * kaonCorr_full[bins_p];

	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_p[i] = new TH3F( Form("hKaonCorrP_%i", i),Form("hKaonCorrP_%i", i), nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
		kaonCorr_m[i] = new TH3F( Form("hKaonCorrM_%i", i),Form("hKaonCorrM_%i", i) , nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
		kaonCorr_full[i] = new TH3F( Form("hKaonCorr_%i", i), "hKaonCorr", nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
	}

	cout<<"BEGIN FITTING\n";
	
	//Loop through bins for fitting
	outFile_fits->cd();
	for( int j = 0; j < nBinsQ2; j++ ){
		for( int k = 0; k < nBinsXb; k++ ){
			for( int l = 0; l < nBinsZ; l++ ){
				for( int m = 0; m < bins_p; m++ ){//momentum bins
					//If not filled, continue
					if(!hBeta[0][j+1][k+1][l+1][m+1]){continue;}
					if(!hBeta[1][j+1][k+1][l+1][m+1]){continue;}
					
					cout<<"Found histogram\n";

					hBeta[0][j+1][k+1][l+1][m+1]->Write();
					hBeta[1][j+1][k+1][l+1][m+1]->Write();
					//fit positive pions

					double fitProb_pip = 0;

					fitPip[j][k][l][m] = new TF1(Form("f_pip_%i_%i_%i_%i", j+1, k+1, l+1, m+1), "gaus", min[m], max[m]);
					fitPip[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum(), .139*.139, .01 ); 
					hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_%i_%i_%i_%i", j+1, k+1, l+1, m+1), "R" );
					fitProb_pip = fitPip[j][k][l][m]->GetProb();	

					//fit negative pions	
					double fitProb_pim = 0;
					fitPim[j][k][l][m] = new TF1(Form("f_pim_%i_%i_%i_%i", j+1, k+1, l+1, m+1), "gaus", min[m], max[m]);
					fitPim[j][k][l][m]->SetParameters( hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum(), .139*.139, .01 ); 
					hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_%i_%i_%i_%i", j+1, k+1, l+1, m+1), "R" );
					fitProb_pim = fitPim[j][k][l][m]->GetProb();	
					
					//Get Pip Counts for correction
					double pip_num = 0;
					double pip_fit_mean = fitPip[j][k][l][m]->GetParameter(1);
					double pip_fit_std = fitPip[j][k][l][m]->GetParameter(2);
					int pip_bin_min = 1;//hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pip_fit_mean - 2*pip_fit_std );
					int pip_bin_max = 1;
					
				 	if( pip_fit_mean < -0.05 || pip_fit_mean > 0.075 || fitProb_pip <= 0 || fitProb_pip >= 1 || hBeta[0][j+1][k+1][l+1][m+1]->Integral() <= 120 ||
							 pip_fit_mean + sigma*pip_fit_std > .2){	
						pip_bin_max= hBeta[0][j+1][k+1][l+1][m+1]->FindBin( 0.1 );
					}
					else{
						pip_bin_max = hBeta[0][j+1][k+1][l+1][m+1]->FindBin( pip_fit_mean + sigma*pip_fit_std );

						fitPip[j][k][l][m]->Write();	
					}
				

					for( int bin = 1; bin <= pip_bin_max; bin++ ){
						pip_num+=hBeta[0][j+1][k+1][l+1][m+1]->GetBinContent(bin);
					}

					double kp_num = hBeta[0][j+1][k+1][l+1][m+1]->Integral(0, hBeta[0][j+1][k+1][l+1][m+1]->GetNbinsX()) - pip_num;
					double pip_err = (1./pow( kp_num+pip_num, 2))*sqrt( pip_num*kp_num*(pip_num+kp_num ) );
					
					//Get Pim Counts for correction
					
					double pim_num = 0;
					double pim_fit_mean = fitPim[j][k][l][m]->GetParameter(1);
					double pim_fit_std = fitPim[j][k][l][m]->GetParameter(2);
					int pim_bin_min = 1;//hBeta[1][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pim_fit_mean - 2*pim_fit_std );
					int pim_bin_max = 1;//hBeta[1][j+1][k+1][l+1][m+1]->GetBin( pim_fit_mean + 2*pim_fit_std );
				 	
					if( pim_fit_mean < 0 || pim_fit_mean > 0.05 || fitProb_pim <= 0 || fitProb_pim >= 1 || hBeta[1][j+1][k+1][l+1][m+1]->Integral() <= 150||
					pim_fit_mean + sigma*pim_fit_std > .2){	
						pim_bin_max= hBeta[1][j+1][k+1][l+1][m+1]->FindBin( 0.1 );
					}
					else{
						pim_bin_max = hBeta[1][j+1][k+1][l+1][m+1]->FindBin( pim_fit_mean + sigma*pim_fit_std );
						fitPim[j][k][l][m]->Write();	
					}

					for( int bin = 1; bin <= pim_bin_max; bin++ ){
						pim_num+=hBeta[1][j+1][k+1][l+1][m+1]->GetBinContent(bin);
					}

					double km_num = hBeta[1][j+1][k+1][l+1][m+1]->Integral(0, hBeta[1][j+1][k+1][l+1][m+1]->GetNbinsX()) - pim_num;
					double pim_err = (1./pow( km_num+pim_num, 2))*sqrt( pim_num*km_num*(pim_num+km_num ) );
				

					if(pip_num/(pip_num+kp_num) >= 0 && pip_num/(pip_num+kp_num) < 1.  ){
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, pip_num/(pip_num+kp_num));
						kaonCorr_p[m]->SetBinError(k+1, j+1, l+1, pip_err);
					}
					else{
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, 0);
						kaonCorr_p[m]->SetBinError(k+1, j+1, l+1, 0);
					}
					if(pim_num/(pim_num+km_num) >= 0  && pim_num/(pim_num+km_num) < 1 ){
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, pim_num/(pim_num+km_num));
						kaonCorr_m[m]->SetBinError(k+1, j+1, l+1, pim_err);
					}
					else{
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, -1);
						kaonCorr_m[m]->SetBinError(k+1, j+1, l+1, 0);
					}
					
					//fitPim[j][k][l][m]->Write();	
					//fitPip[j][k][l][m]->Write();	
			
				}
			}
		}
	}	
	outFile_fits->Close();
	outFile->cd();
	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_full[i] = (TH3F*)kaonCorr_p[i]->Clone(Form("hKaonCorr_%i", i)); 
		kaonCorr_full[i]->Divide(kaonCorr_m[i]);
	      	kaonCorr_p[i]->Write();
		kaonCorr_m[i]->Write();
		kaonCorr_full[i]->Write();
		kaonCorr_full[i]->Print();
	}	

	outFile->Close();
}

