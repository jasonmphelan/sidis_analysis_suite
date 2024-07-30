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

       	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
       	TFile * inFile = new TFile((TString) HIST_PATH + "/" + in_name + ".root");


	
	cout<<"GETTING HISTS\n";

	TH1F * hBeta[2][bins_Q2+1][bins_xB+1][bins_Z+1][bins_p+1];

	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){
			for( int l = 0; l <= bins_Z; l++ ){
				for( int m = 0; m <= bins_p; m++ ){//momentum bins
					hBeta[0][j][k][l][m] = (TH1F*)inFile->Get((TString)"hBeta_rich_pip"+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m)); 
					hBeta[1][j][k][l][m] = (TH1F*)inFile->Get((TString)"hBeta_rich_pim"+Form("_Q2_%i_xB_%i_Z_%i_p_%i", j, k, l, m)); 
					
				}
			}
		}
	}

	TF1 *fitPim[bins_Q2][bins_xB][bins_Z][bins_p];
	TF1 *fitPip[bins_Q2][bins_xB][bins_Z][bins_p];
	TF1 *fitPim_k[bins_Q2][bins_xB][bins_Z][bins_p];
	TF1 *fitPip_k[bins_Q2][bins_xB][bins_Z][bins_p];
	
	double min[4] = { .995, .995, .995, .995 };
	double max[4] = { 1.05, 1.05, 1.05,  1.05 };
	double min_k[4] = { .98, .98, .98, .98 };

	TH3F * kaonCorr_p[bins_p];
	TH3F * kaonCorr_m[bins_p];
	TH3F * kaonCorr_full[bins_p];

	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_p[i] = new TH3F( Form("hKaonCorrP_%i", i), "hKaonCorrP", bins_xB, xB_min, xB_max, bins_Q2, Q2_min, Q2_max, bins_Z, Z_min, Z_max); 
		kaonCorr_m[i] = new TH3F( Form("hKaonCorrM_%i", i), "hKaonCorrM", bins_xB, xB_min, xB_max, bins_Q2, Q2_min, Q2_max, bins_Z, Z_min, Z_max); 
		kaonCorr_full[i] = new TH3F( Form("hKaonCorr_%i", i), "hKaonCorr", bins_xB, xB_min, xB_max, bins_Q2, Q2_min, Q2_max, bins_Z, Z_min, Z_max); 
	}

	cout<<"BEGIN FITTING\n";
	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();
	for( int j = 0; j < bins_Q2; j++ ){
		for( int k = 0; k < bins_xB; k++ ){
			for( int l = 0; l < bins_Z; l++ ){
				for( int m = 0; m < bins_p; m++ ){//momentum bins

					if(!hBeta[0][j+1][k+1][l+1][m+1]){continue;}
					if(!hBeta[1][j+1][k+1][l+1][m+1]){continue;}
					
					//fit positive pions

					fitPip[j][k][l][m] = new TF1(Form("f_pip_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
					fitPip[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[0][j+1][k+1][l+1][m+1]->GetMean(), hBeta[0][j+1][k+1][l+1][m+1]->GetRMS() ); 
					hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_%i_%i_%i_%i", j, k, l, m), "R" );
					
					//fit negative pions	
					fitPim[j][k][l][m] = new TF1(Form("f_pim_%i_%i_%i_%i", j, k, l, m), "gaus", min[m], max[m]);
					fitPim[j][k][l][m]->SetParameters( hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum(), hBeta[1][j+1][k+1][l+1][m+1]->GetMean(), hBeta[1][j+1][k+1][l+1][m+1]->GetRMS() ); 
					hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_%i_%i_%i_%i", j, k, l, m), "R" );
					//fit positive kaons

					fitPip_k[j][k][l][m] = new TF1(Form("f_pip_k_%i_%i_%i_%i", j, k, l, m), "gaus", min_k[m], min[m]);
					fitPip_k[j][k][l][m]->SetParameters( hBeta[0][j+1][k+1][l+1][m+1]->GetMaximum()/5., min[m]-min_k[m], hBeta[0][j+1][k+1][l+1][m+1]->GetRMS()/2. ); 
					hBeta[0][j+1][k+1][l+1][m+1]->Fit( Form("f_pip_k_%i_%i_%i_%i", j, k, l, m),"R" );
					
					//fit positive kaons

					fitPim_k[j][k][l][m] = new TF1(Form("f_pim_k_%i_%i_%i_%i", j, k, l, m), "gaus", min_k[m], min[m]);
					fitPim_k[j][k][l][m]->SetParameters( hBeta[1][j+1][k+1][l+1][m+1]->GetMaximum()/5., min[m]-min_k[m], hBeta[1][j+1][k+1][l+1][m+1]->GetRMS()/2. ); 
					hBeta[1][j+1][k+1][l+1][m+1]->Fit( Form("f_pim_k_%i_%i_%i_%i", j, k, l, m), "R" );
					
					//fitPip[j][k][l][m]->SetRange(.9, 1.1); 
					//fitPim[j][k][l][m]->SetRange(.9, 1.1); 
					//fitPip_k[j][k][l][m]->SetRange(.9, 1.1); 
					//fitPim_k[j][k][l][m]->SetRange(.9, 1.1); 

					double pip_num = fitPip[j][k][l][m]->Integral(.9, 1.1);
					double pip_den = pip_num + fitPip_k[j][k][l][m]->Integral(.9, 1.1);
					
					double pim_num = fitPim[j][k][l][m]->Integral(.9, 1.1);
					double pim_den = pim_num + fitPim_k[j][k][l][m]->Integral(.9, 1.1);
					
					if(pip_num/pip_den > .5){
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, pip_num/pip_den);
					}
					else{
						kaonCorr_p[m]->SetBinContent(k+1, j+1, l+1, 0);
					}
					if(pim_num/pim_den > .5){
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, pim_num/pim_den);
					}
					else{
						kaonCorr_m[m]->SetBinContent(k+1, j+1, l+1, -1);
					}
					/*
					double pip_corr = (double) fitPip[j][k][l][m]->Integral(.9, 1.1) / ( hBeta[0][j+1][k+1][l+1][m+1]->GetBinWidth(1)*(double) (hBeta[0][j+1][k+1][l+1][m+1]->Integral() ) );
					double pim_corr = (double) fitPim[j][k][l][m]->Integral(.9, 1.1) / ( hBeta[1][j+1][k+1][l+1][m+1]->GetBinWidth(1)*(double) (hBeta[1][j+1][k+1][l+1][m+1]->Integral() ) );
					if(pip_corr > 1 || pim_corr > 1){
						kaonCorr_p[m]->SetBinContent( k+1, j+1, l+1, 0); 

					//fitPim[j][k][l][m]->Write();			
						kaonCorr_m[m]->SetBinContent( k+1, j+1, l+1, -1);
					}
					else{
						kaonCorr_p[m]->SetBinContent( k+1, j+1, l+1, pip_corr); 

					//fitPim[j][k][l][m]->Write();			
						kaonCorr_m[m]->SetBinContent( k+1, j+1, l+1, pim_corr);
					}
					*/		
					
					cout<<"Pip Value : "<<kaonCorr_p[m]->GetBinContent(k+1, j+1, l+1)<<std::endl;	
					cout<<"Pim Value : "<<kaonCorr_m[m]->GetBinContent(k+1, j+1, l+1)<<std::endl;	
					//hBeta[0][p]->SetMarkerStyle(kCircle);
					hBeta[0][j+1][k+1][l+1][m+1]->Draw("");
					fitPip[j][k][l][m]->Draw("SAME");	
					fitPip_k[j][k][l][m]->Draw("SAME");	
					canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");
					canvas.Clear();
						
						
					//hBeta[1][p]->SetMarkerStyle(kCircle);
					hBeta[1][j+1][k+1][l+1][m+1]->Draw("");
					fitPim[j][k][l][m]->Draw("SAME");	
					fitPim_k[j][k][l][m]->Draw("SAME");	
					canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");
					canvas.Clear();
			
				}
			}
		}
	}	
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
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


double fitF( double *p, double *par){
	return par[0] + par[1]*p[0] + par[2]*p[0] + par[3]*p[0];
}
