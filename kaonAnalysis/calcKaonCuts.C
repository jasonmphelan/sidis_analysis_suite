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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

const int NpBins = 4;

const double p_min = 1.25;
const double p_max = 4.75;
//const double Q2_width = (Q2_max - Q2_min)/( (double) nQ2Bins);

double fitF( double *p, double *par);

void calcKaonCuts(){


	TFile * outFile = new TFile("kaonCutFunc.root", "RECREATE");
	TFile inFile("../histMakers/kinematics/kinematics_data_10.2_beta.root");
	//Get Cut Value (cut out 10% percentile)
	
	outFile->cd();
	
	cout<<"GETTING HISTS\n";

	TH1F * hBeta[2][4]; 	

	for( int p = 1; p <= 4; p++ ){
		hBeta[0][p-1] = (TH1F*)inFile.Get(Form("hBeta_rich_pip_%i_0", p)); 
		hBeta[1][p-1] = (TH1F*)inFile.Get(Form("hBeta_rich_pim_%i_0", p)); 
	}	

	TF1 *fitPim[4];
	TF1 *fitPip[4];
	TF1 *fitTheta[4];
	
	double minPip[4] = { .992, .998, .998, .999 };
	double maxPip[4] = { 1.01, 1.01, 1.01,  1.01 };
	
	double minPim[4] = { .992, .994, .998, .999 };
	double maxPim[4] = { 1.01, 1.01, 1.01,  1.01 };

	cout<<"BEGIN FITTING\n";
	for( int p = 0; p < 4; p++ ){


		//fitPip[p] = new TF1(Form("f_pip_%i", p), "[0]+[1]*x + [2]*x*x + [3]*x*x*x", 1.25, 4.75);
		fitPip[p] = new TF1(Form("f_pip_%i", p), "gaus", minPip[p], maxPip[p]);
		fitPip[p]->SetParameters( hBeta[0][p]->GetMaximum(), hBeta[0][p]->GetMean(), hBeta[0][p]->GetRMS() ); 

		hBeta[0][p]->Fit( Form("f_pip_%i", p) );
		
		fitPim[p] = new TF1(Form("f_pim_%i", p), "gaus", minPim[p], maxPim[p]);
		fitPim[p]->SetParameters( hBeta[1][p]->GetMaximum(), hBeta[1][p]->GetMean(), hBeta[1][p]->GetRMS() ); 

		hBeta[1][p]->Fit( Form("f_pim_%i", p) );

		fitPip[p]->Write();			
		fitPim[p]->Write();			
		
		TCanvas * c1 = new TCanvas( Form("cPip_%i", p) );
		//hBeta[0][p]->SetMarkerStyle(kCircle);
		hBeta[0][p]->Draw("");
		fitPip[p]->Draw("SAME");	
		
		
		TCanvas * c2 = new TCanvas( Form("cPim_%i", p) );
		//hBeta[1][p]->SetMarkerStyle(kCircle);
		hBeta[1][p]->Draw("");
		fitPim[p]->Draw("SAME");	
		
		c1->Write();
		c2->Write();	
	}
	
	outFile->Close();
}


double fitF( double *p, double *par){
	return par[0] + par[1]*p[0] + par[2]*p[0] + par[3]*p[0];
}
