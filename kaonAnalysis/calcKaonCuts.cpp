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

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File (no extension)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

       	TFile * outFile = new TFile((TString) CORR_PATH + "/correctionFiles/" + out_name + ".root", "RECREATE");
       	TFile * inFile = new TFile((TString) HIST_PATH + "/" + in_name + ".root");
	
	int nBinsQ2 = bins_Q2/2;
	int nBinsXb = bins_xB/2;
	int nBinsZ = 2*bins_Z;


	
	cout<<"GETTING HISTS\n";

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

	TF1 *fitPim[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	TF1 *fitPip[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	TF1 *fitPim_k[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	TF1 *fitPip_k[nBinsQ2][nBinsXb][nBinsZ][bins_p];
	
	double min[4] = { .995, .995, .995, .995 };
	double max[4] = { 1.005, 1.005, 1.005,  1.005 };
	double min_k[4] = { .985, .985, .985, .985 };

	TH3F * kaonCorr_p[bins_p];
	TH3F * kaonCorr_m[bins_p];
	TH3F * kaonCorr_full[bins_p];

	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_p[i] = new TH3F( Form("hKaonCorrP_%i", i),Form("hKaonCorrP_%i", i), nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
		kaonCorr_m[i] = new TH3F( Form("hKaonCorrM_%i", i),Form("hKaonCorrM_%i", i) , nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
		kaonCorr_full[i] = new TH3F( Form("hKaonCorr_%i", i), "hKaonCorr", nBinsXb, xB_min, xB_max, nBinsQ2, Q2_min, Q2_max, nBinsZ, Z_min, Z_max); 
	}

	cout<<"BEGIN FITTING\n";
	TCanvas canvasP("canvasP");
	canvasP.Print((TString) HIST_PATH + "/" + out_name + "_p.pdf[");
	canvasP.Clear();
	TCanvas canvasM("canvaMP");
	canvasM.Print((TString) HIST_PATH + "/" + out_name + "_m.pdf[");
	canvasM.Clear();
	for( int j = 0; j < nBinsQ2; j++ ){
		for( int k = 0; k < nBinsXb; k++ ){
			for( int l = 0; l < nBinsZ; l++ ){
				for( int m = 0; m < bins_p; m++ ){//momentum bins

					if(!hBeta[0][j+1][k+1][l+1][m+1]){continue;}
					if(!hBeta[1][j+1][k+1][l+1][m+1]){continue;}
					
					//fit positive pions

					double fitProb = 0;

					fitPip[j][k][l][m] = new TF1(Form("f_pip_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
					fitPip[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[0][j+1][k+1][l+1][m+1]->GetMean(), hBeta[0][j+1][k+1][l+1][m+1]->GetRMS() ); 
					hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_%i_%i_%i_%i", j, k, l, m), "R" );
					fitProb = fitPip[j][k][l][m]->GetProb();	
				
					//if bad fit, try lorentzian	
					if( fitPip[j][k][l][m]->GetProb() < .75 ){
						fitPip[j][k][l][m] = new TF1(Form("f_pip_%i_%i_%i_%i", j, k, l, m), "gaus",
										min[m]+.0005, max[m]);
										//[0] * (1./TMath::Pi()) * (1./[2]) / ( pow ( (x - [1]), 2 ) + pow( [2], 2 ) ), 
						
						fitPip[j][k][l][m]->SetParNames("Maximum", "Mean", "Width");
						fitPip[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[0][j+1][k+1][l+1][m+1]->GetMean(), hBeta[0][j+1][k+1][l+1][m+1]->GetRMS() ); 
						hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_%i_%i_%i_%i", j, k, l, m), "R" );
					}
					//if no luck, go back to gauss
					if( fitProb > fitPip[j][k][l][m]->GetProb() ){
						fitPip[j][k][l][m] = new TF1(Form("f_pip_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
						fitPip[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[0][j+1][k+1][l+1][m+1]->GetMean(), hBeta[0][j+1][k+1][l+1][m+1]->GetRMS() ); 
						hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_%i_%i_%i_%i", j, k, l, m), "R" );
					}
					
						
					//fit negative pions	
					fitPim[j][k][l][m] = new TF1(Form("f_pim_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
					fitPim[j][k][l][m]->SetParameters( hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[1][j+1][k+1][l+1][m+1]->GetMean(), hBeta[1][j+1][k+1][l+1][m+1]->GetRMS() ); 
					hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_%i_%i_%i_%i", j, k, l, m), "R" );
					
					fitProb = fitPim[j][k][l][m]->GetProb();	
					
					if( fitPim[j][k][l][m]->GetProb() < .75 ){
						fitPim[j][k][l][m] = new TF1(Form("f_pim_%i_%i_%i_%i", j, k, l, m), 
										"gaus",
										min[m] + .0005, max[m]);
										//[0] * (1./TMath::Pi()) * (1./[2]) / ( pow ( (x - [1]), 2 ) + pow( [2], 2 ) ), 
						
						fitPim[j][k][l][m]->SetParNames("Maximum", "Mean", "Width");
						fitPim[j][k][l][m]->SetParameters( hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[1][j+1][k+1][l+1][m+1]->GetMean(), hBeta[1][j+1][k+1][l+1][m+1]->GetRMS() ); 
						hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_%i_%i_%i_%i", j, k, l, m), "R" );
					}
					if( fitProb > fitPim[j][k][l][m]->GetProb() ){
						fitPim[j][k][l][m] = new TF1(Form("f_pim_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
						fitPim[j][k][l][m]->SetParameters( 
									hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum(), 
									hBeta[1][j+1][k+1][l+1][m+1]->GetMean(), 
									hBeta[1][j+1][k+1][l+1][m+1]->GetRMS() ); 
						
						hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_%i_%i_%i_%i", j, k, l, m), "R" );
					
					}	
					
					double pip_num = 0;
					double pip_fit_mean = fitPip[j][k][l][m]->GetParameter(1);
					double pip_fit_std = fitPip[j][k][l][m]->GetParameter(2);
					int pip_bin_min = hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pip_fit_mean - 2*pip_fit_std );
					//int pip_bin_max = hBeta[0][j+1][k+1][l+1][m+1]->GetBin( pip_fit_mean + 2*pip_fit_std );
					int pip_bin_max = hBeta[0][j+1][k+1][l+1][m+1]->GetNbinsX( );
					
					for( int bin = pip_bin_min; bin <= pip_bin_max; bin++ ){
						pip_num+=hBeta[0][j+1][k+1][l+1][m+1]->GetBinContent(bin);
					}

					double kp_num = hBeta[0][j+1][k+1][l+1][m+1]->Integral() - pip_num;
					double pip_err = (1./pow( kp_num+pip_num, 2))*sqrt( pip_num*kp_num*(pip_num+kp_num ) );
					
					double pim_num = 0;
					double pim_fit_mean = fitPim[j][k][l][m]->GetParameter(1);
					double pim_fit_std = fitPim[j][k][l][m]->GetParameter(2);
					int pim_bin_min = hBeta[1][j+1][k+1][l+1][m+1]->GetXaxis()->FindBin( pim_fit_mean - 2*pim_fit_std );
					int pim_bin_max = hBeta[1][j+1][k+1][l+1][m+1]->GetNbinsX( );
					
					for( int bin = pim_bin_min; bin <= pim_bin_max; bin++ ){
						pim_num+=hBeta[1][j+1][k+1][l+1][m+1]->GetBinContent(bin);
					}

					double km_num = hBeta[1][j+1][k+1][l+1][m+1]->Integral() - pim_num;
					double pim_err = (1./pow( km_num+pim_num, 2))*sqrt( pim_num*km_num*(pim_num+km_num ) );
					
					cout<<"MIN BIN : "<<pim_bin_min<<std::endl;
					cout<<"KAONS : "<<km_num<<std::endl;
					cout<<"FIT PROBABILITY : "<<fitPim[j][k][l][m]->GetProb()<<"\n";	
					cout<<"PIM CORRECTION : "<<pim_num/(pim_num+km_num)<<std::endl;


					if(pip_num/(pip_num+kp_num) > .25 && pip_num/(pip_num+kp_num) < 1. ){
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, pip_num/(pip_num+kp_num));
						kaonCorr_p[m]->SetBinError(k+1, j+1, l+1, pip_err);
					}
					else{
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, 0);
						kaonCorr_p[m]->SetBinError(k+1, j+1, l+1, 0);
					}
					if(pim_num/(pim_num+km_num) > .25  && pim_num/(pim_num+km_num) < 1 ){
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, pim_num/(pim_num+km_num));
						kaonCorr_m[m]->SetBinError(k+1, j+1, l+1, pim_err);
					}
					else{
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, -1);
						kaonCorr_m[m]->SetBinError(k+1, j+1, l+1, 0);
					}
					
					if( j == 0 && k == 1 && l == 4 && m == 2 ){
						canvasP.cd();
						hBeta[0][j+1][k+1][l+1][m+1]->SetStats(0);;
						hBeta[0][j+1][k+1][l+1][m+1]->SetTitle("");
						hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->SetTitle("#beta");
						hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->SetTitleFont(43);
        					hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->SetTitleSize(20);
        					hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->SetRangeUser(.985, 1.005);
        					hBeta[0][j+1][k+1][l+1][m+1]->GetXaxis()->SetTitleOffset(0.9);
						hBeta[0][j+1][k+1][l+1][m+1]->GetYaxis()->SetTitleFont(43);
        					hBeta[0][j+1][k+1][l+1][m+1]->GetYaxis()->SetTitleSize(20);
        					hBeta[0][j+1][k+1][l+1][m+1]->GetYaxis()->SetTitleOffset(1.25);
						hBeta[0][j+1][k+1][l+1][m+1]->Draw("");

						fitPip[j][k][l][m]->Draw("SAME");	
						canvasP.Print((TString) HIST_PATH + "/" + out_name + "_p.pdf");
						canvasP.Clear();
					}
						
					//hBeta[1][j+1][k+1][l+1][m+1]->Draw("");
					//canvasM.cd();
					//fitPim[j][k][l][m]->Draw("SAME");	
					//canvasM.Print((TString) HIST_PATH + "/" + out_name + ".pdf_m");
					//canvasM.Clear();
			
				}
			}
		}
	}	
	canvasP.Print((TString) HIST_PATH + "/" + out_name + "_p.pdf]");
	canvasM.Print((TString) HIST_PATH + "/" + out_name + "_m.pdf]");
	outFile->cd();
	
	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_full[i] = (TH3F*)kaonCorr_p[i]->Clone(Form("hKaonCorr_%i", i)); 
		kaonCorr_full[i]->Divide(kaonCorr_m[i]);
	      	kaonCorr_p[i]->Write();
		kaonCorr_m[i]->Write();
		kaonCorr_full[i]->Write();
	}	

	outFile->Close();
}

