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
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

       	TFile * outFile = new TFile( out_name, "RECREATE");
       	TFile * inFile = new TFile(in_name);
	
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

	for( int p = 0; p <  bins_p; p++ ){
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
						p_fit_min = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) - .025/2.;
					}
					if( p_fit_max == 1 && p_fit_min > 0 &&  kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						p_fit_max = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) - .025/2.;
					}
					if( kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						kaonCorr_m_1d[p][q][x]->SetBinContent( z+1, kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ));
						kaonCorr_m_1d[p][q][x]->SetBinError( z+1, kaonCorr_m[p]->GetBinError( x+1, q+1, z+1 ));
					}
					if( m_fit_min == 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						m_fit_min = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) - .025/2.;
					}
					if( m_fit_max == 1 && m_fit_min > 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						m_fit_max = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) - .025/2.;
				
					}
				}

				kaonCorr_p_1d[p][q][x]->Write();
				kaonCorr_m_1d[p][q][x]->Write();

				fitPip[p][q][x] = new TF1( Form("fitPip_%i_%i_%i", p, q, x), "[0]+ [1]*x", p_fit_min, p_fit_max );
				fitPip[p][q][x]->SetParameters( .5, .5);
				
				fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*x", m_fit_min, m_fit_max );
				fitPim[p][q][x]->SetParameters( 0.5, 0.5 );
			
				if( kaonCorr_p_1d[p][q][x]->GetEntries() <= 1 ){
					fitPip[p][q][x]->SetParameters(  0, 0, 0 );
				}
					
				else{
					kaonCorr_p_1d[p][q][x]->GetXaxis()->SetRangeUser( p_fit_min, p_fit_max );
					kaonCorr_p_1d[p][q][x]->Smooth(1, "R");
					kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				}
				
				if( kaonCorr_m_1d[p][q][x]->GetEntries() <= 1 ){
					fitPim[p][q][x]->SetParameters( 0, 0, 0 );
				}
					
				else{
					kaonCorr_m_1d[p][q][x]->GetXaxis()->SetRangeUser( m_fit_min, m_fit_max );
					kaonCorr_m_1d[p][q][x]->Smooth(1, "R");
					kaonCorr_m_1d[p][q][x]->Fit( Form("fitPim_%i_%i_%i", p, q, x) );
				}
				
				
				kaonCorr_p_1d[p][q][x]->Write(Form("kaonCorr_p_1d_%i_%i_%i_smooth", p, q, x));
				kaonCorr_m_1d[p][q][x]->Write(Form("kaonCorr_m_1d_%i_%i_%i_smooth", p, q, x));
				
				fitPim[p][q][x]->Write();		
				fitPip[p][q][x]->Write();		
			}

		}
	}
	


	outFile->Close();
}

