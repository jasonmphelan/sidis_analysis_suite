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
       	TFile * inFile = new TFile((TString) CORR_PATH + "/correctionFiles/" + in_name + ".root");
	
	int nBinsQ2 = bins_Q2;
	int nBinsXb = bins_xB;
	int nBinsZ = 2*bins_Z;

	
	cout<<"GETTING HISTS\n";

	TH3F * kaonCorr_p[bins_p];
	TH3F * kaonCorr_m[bins_p];

	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_p[i] = (TH3F *)inFile->Get( Form("hKaonCorrP_%i", i) );
		kaonCorr_m[i] = (TH3F *)inFile->Get( Form("hKaonCorrM_%i", i) );
	}
	
	//Project to 1D histograms
	TH1F * kaonCorr_p_1d[bins_p][nBinsQ2][nBinsXb];
	TH1F * kaonCorr_m_1d[bins_p][nBinsQ2][nBinsXb];
	TF1 *fitPim[nBinsQ2][nBinsXb][bins_p];
	TF1 *fitPip[nBinsQ2][nBinsXb][bins_p];

	outFile->cd();
	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int p = 2; p <  bins_p; p++ ){
		for( int x = 0; x < nBinsXb; x++ ){
			for( int q = 0; q <  nBinsQ2; q++ ){
				kaonCorr_p_1d[p][q][x] = new TH1F( Form("kaonCorr_p_1d_%i_%i_%i", p, q, x), ";z;w+", nBinsZ, Z_min, Z_max );
				kaonCorr_m_1d[p][q][x] = new TH1F( Form("kaonCorr_m_1d_%i_%i_%i", p, q, x), ";z;w-", nBinsZ, Z_min, Z_max );
			
				double p_fit_min = 0;
				double p_fit_max = 1;
				double m_fit_min = 0;
				double m_fit_max = 1;

				for( int z = 0; z < nBinsZ; z++ ){
					if( kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){	
						kaonCorr_p_1d[p][q][x]->SetBinContent( z+1, kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ));
						kaonCorr_p_1d[p][q][x]->SetBinError( z+1, kaonCorr_p[p]->GetBinError( x+1, q+1, z+1 ));
					}

					if( p_fit_min == 0 &&  kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						p_fit_min = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) - .025;
					}
					if( p_fit_max == 1 && p_fit_min > 0 &&  kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						p_fit_max = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) + .025;
					}

					kaonCorr_m_1d[p][q][x]->SetBinContent( z+1, kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_m_1d[p][q][x]->SetBinError( z+1, kaonCorr_m[p]->GetBinError( x+1, q+1, z+1 ));
					
					if( m_fit_min == 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						m_fit_min = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) - .025;
					}
					if( m_fit_max == 1 && m_fit_min > 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						m_fit_max = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) + .025;
					}
				}

				fitPip[p][q][x] = new TF1( Form("fitPip_%i_%i_%i", p, q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", p_fit_min, p_fit_max );
				//fitPip[p][q][x] = new TF1( Form("fitPip_%i_%i_%i", p, q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
				fitPip[p][q][x]->SetParameters( 1, 1, 1, 1);
				
				fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", m_fit_min, m_fit_max );
				//fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
				fitPim[p][q][x]->SetParameters( 1, 1, 1, 1 );

				kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				kaonCorr_m_1d[p][q][x]->Fit( Form("fitPim_%i_%i_%i", p, q, x) );
				
				/*		
				if( fitPip[p][q][x]->Eval(1) > 0 ){
					p_fit_max = 1;
					kaonCorr_p_1d[p][q][x]->SetBinContent(nBinsZ, .5);			
					kaonCorr_p_1d[p][q][x]->SetBinError(nBinsZ, .5);			
					kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				}
				if( fitPip[p][q][x]->Eval(1) < 0 ){
					p_fit_max = 1;
					kaonCorr_p_1d[p][q][x]->SetBinContent(nBinsZ, .5);			
					kaonCorr_p_1d[p][q][x]->SetBinError(nBinsZ, 0.5);			
				}
				if( fitPip[p][q][x]->Eval(0) > 0 ){
					p_fit_min = 0.3;
					fitPip[p][q][x]->SetRange(.3, p_fit_max);
					kaonCorr_p_1d[p][q][x]->SetBinContent(1, .5);			
					kaonCorr_p_1d[p][q][x]->SetBinError(1, .5);			
				}
				if( fitPip[p][q][x]->Eval(0) < 0 ){
					p_fit_min = 0.3;
					kaonCorr_p_1d[p][q][x]->SetBinContent(1, 0.5);			
					kaonCorr_p_1d[p][q][x]->SetBinError(1, 0.5);			
				}
				if( fitPim[p][q][x]->Eval(1) > 0 ){
					m_fit_max = 1;
					kaonCorr_m_1d[p][q][x]->SetBinContent(nBinsZ, .5);			
					kaonCorr_m_1d[p][q][x]->SetBinError(nBinsZ, .5);			
				}
				if( fitPim[p][q][x]->Eval(1) < 0 ){
					m_fit_max = 1;
					kaonCorr_m_1d[p][q][x]->SetBinContent(nBinsZ, 0.5);			
					kaonCorr_m_1d[p][q][x]->SetBinError(nBinsZ, 0.5);			
				}
				if( fitPim[p][q][x]->Eval(0) > 0 ){
					m_fit_min = .3;
					kaonCorr_m_1d[p][q][x]->SetBinContent(1, .5);			
					kaonCorr_m_1d[p][q][x]->SetBinError(1, .5);			
				}
				if( fitPim[p][q][x]->Eval(0) < 0 ){
					m_fit_min = .3;
					kaonCorr_m_1d[p][q][x]->SetBinContent(1, 0.5);			
					kaonCorr_m_1d[p][q][x]->SetBinError(1, 0.5);			
				}
				fitPip[p][q][x]->SetRange(p_fit_min, p_fit_max);
				fitPim[p][q][x]->SetRange(m_fit_min, m_fit_max);
					
				kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				kaonCorr_m_1d[p][q][x]->Fit( Form("fitPim_%i_%i_%i", p, q, x) );
				*/

				kaonCorr_p_1d[p][q][x]->SetLineColor(kRed);
				kaonCorr_p_1d[p][q][x]->GetYaxis()->SetRangeUser(0, 1);
				kaonCorr_p_1d[p][q][x]->Draw();
				//fitPip[p][q][x]->SetLineColor(kBlue);
				fitPip[p][q][x]->Draw("SAME");
				//kaonCorr_m_1d[p][q][x]->SetLineColor(kRed);
				kaonCorr_m_1d[p][q][x]->Draw("SAME");
				fitPim[p][q][x]->SetLineColor(kMagenta);
				
				fitPim[p][q][x]->Draw("SAME");
				canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
				canvas.Clear();
			
			}

		}
	}
	


	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
	outFile->Close();
}
