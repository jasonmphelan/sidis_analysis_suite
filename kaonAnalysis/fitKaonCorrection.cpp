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


	
	cout<<"GETTING HISTS\n";

	TH3F * kaonCorr_p[bins_p];
	TH3F * kaonCorr_m[bins_p];

	for( int i = 0; i < bins_p; i++ ){
		kaonCorr_p[i] = (TH3F *)inFile->Get( Form("hKaonCorrP_%i", i) );
		kaonCorr_m[i] = (TH3F *)inFile->Get( Form("hKaonCorrM_%i", i) );
	}
	
	//Project to 1D histograms
	TH1F * kaonCorr_p_1d[bins_p][bins_Q2][bins_xB];
	TH1F * kaonCorr_m_1d[bins_p][bins_Q2][bins_xB];
	TF1 *fitPim[bins_Q2][bins_xB][bins_p];
	TF1 *fitPip[bins_Q2][bins_xB][bins_p];
	
	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int p = 0; p <  bins_p; p++ ){
		for( int q = 0; q <  bins_Q2; q++ ){
			for( int x = 0; x < bins_xB; x++ ){
				kaonCorr_p_1d[p][q][x] = new TH1F( Form("kaonCorr_p_1d_%i_%i_%i", p, q, x), ";z;w+", bins_Z, Z_min, Z_max );
				kaonCorr_m_1d[p][q][x] = new TH1F( Form("kaonCorr_m_1d_%i_%i_%i", p, q, x), ";z;w-", bins_Z, Z_min, Z_max );
			
				double p_fit_min = 0;
				double p_fit_max = 1;
				double m_fit_min = 0;
				double m_fit_max = 1;

				for( int z = 0; z < bins_Z; z++ ){
					kaonCorr_p_1d[p][q][x]->SetBinContent( z+1, kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_p_1d[p][q][x]->SetBinError( z+1, kaonCorr_p[p]->GetBinError( x+1, q+1, z+1 ));
					
					if( p_fit_min == 0 &&  kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						p_fit_min = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) - .025;
					}
					if( p_fit_max == 1 && p_fit_min > 0 &&  kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						p_fit_max = kaonCorr_p_1d[p][q][x]->GetBinCenter(z+1) + .025;
					}

					std::cout<<"min : "<<p_fit_min<<" && max : "<<p_fit_max<<std::endl;
						
					kaonCorr_m_1d[p][q][x]->SetBinContent( z+1, kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_m_1d[p][q][x]->SetBinError( z+1, kaonCorr_m[p]->GetBinError( x+1, q+1, z+1 ));
					
					if( m_fit_min == 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) > 0 ){
						m_fit_min = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) - .025;
					}
					if( m_fit_max == 1 && m_fit_min > 0 &&  kaonCorr_m[p]->GetBinContent( x+1, q+1, z+1 ) <= 0){
						m_fit_max = kaonCorr_m_1d[p][q][x]->GetBinCenter(z+1) + .025;
					}
				}

				//fitPip[p][q][x] = new TF1( Form("fitPip_%i_%i_%i", p, q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", Z_min, Z_max );
				fitPip[p][q][x] = new TF1( Form("fitPip_%i_%i_%i", p, q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
				fitPip[p][q][x]->SetParameters( 1, 1, 1, 1 );
				
				//fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", Z_min, Z_max );
				fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
				fitPim[p][q][x]->SetParameters( 1, 1, 1, 1 );

				kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				kaonCorr_m_1d[p][q][x]->Fit( Form("fitPim_%i_%i_%i", p, q, x) );
						
				kaonCorr_p_1d[p][q][x]->Draw();
				fitPip[p][q][x]->Draw("SAME");
				canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
				canvas.Clear();

			}
		}
	}
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
	
}

