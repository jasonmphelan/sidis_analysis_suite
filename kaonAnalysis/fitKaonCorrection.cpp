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
	int nBinsZ = bins_Z*2;

	
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

	//Fit parameters as a function of Q2 for fixed x
	TH1F * hFitPipQ[bins_p][nBinsXb][4];
	TH1F * hFitPimQ[bins_p][nBinsXb][4];
	TF1 * fitPipQ[bins_p][nBinsXb][4];
	TF1 * fitPimQ[bins_p][nBinsXb][4];

	//Fit parameters as a function of x
	TH1F * hFitPipX[bins_p][4][4];
	TH1F * hFitPimX[bins_p][4][4];
	
	TF1 * fitPipX[bins_p][4][4];
	TF1 * fitPimX[bins_p][4][4];

	for( int p = 0; p < bins_p; p++ ){
		for( int q = 0; q < bins_xB; q++ ){
			for( int c = 0; c < 4; c++ ){
				hFitPipQ[p][q][c] = new TH1F(Form("hFitPipQ_%i_%i_%i", p, q, c ), "", bins_Q2, Q2_min, Q2_max);
				hFitPimQ[p][q][c] = new TH1F(Form("hFitPimQ_%i_%i_%i", p, q, c ), "", bins_Q2, Q2_min, Q2_max);

			}
		}
		for( int c = 0; c < 4; c++ ){
			for( int k = 0; k < 4; k++ ){
				hFitPipX[p][c][k] = new TH1F(Form("hFitPipX_%i_%i_%i", p, c, k ), "", bins_xB, xB_min, xB_max);
				hFitPimX[p][c][k] = new TH1F(Form("hFitPimX_%i_%i_%i", p, c, k ), "", bins_xB, xB_min, xB_max);
			}
		}
	}
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
					kaonCorr_p_1d[p][q][x]->SetBinContent( z+1, kaonCorr_p[p]->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_p_1d[p][q][x]->SetBinError( z+1, kaonCorr_p[p]->GetBinError( x+1, q+1, z+1 ));
					
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
				fitPip[p][q][x]->SetParameters( 1, 1, 1, 1 );
				
				fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", m_fit_min, m_fit_max );
				//fitPim[p][q][x] = new TF1( Form("fitPim_%i_%i_%i", p, q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
				fitPim[p][q][x]->SetParameters( 1, 1, 1, 1 );

				kaonCorr_p_1d[p][q][x]->Fit( Form("fitPip_%i_%i_%i", p, q, x) );
				kaonCorr_m_1d[p][q][x]->Fit( Form("fitPim_%i_%i_%i", p, q, x) );
						
				kaonCorr_p_1d[p][q][x]->SetLineColor(kAzure);
				kaonCorr_p_1d[p][q][x]->Draw();
				fitPip[p][q][x]->SetLineColor(kBlue);
				fitPip[p][q][x]->Draw("SAME");
				kaonCorr_m_1d[p][q][x]->SetLineColor(kRed);
				kaonCorr_m_1d[p][q][x]->Draw("SAME");
				fitPim[p][q][x]->SetLineColor(kMagenta);
				
				fitPim[p][q][x]->Draw("SAME");
				canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
				canvas.Clear();
			
				for( int c = 0; c < 4; c++ ){
					if( fitPip[p][q][x]->GetParError(c) != 0 ){
						hFitPipQ[p][x][c]->SetBinContent( q+1, fitPip[p][q][x]->GetParameter(c));
						hFitPipQ[p][x][c]->SetBinError( q+1, fitPip[p][q][x]->GetParError(c));
					}
					if( fitPim[p][q][x]->GetParError(c) != 0 ){
						hFitPimQ[p][x][c]->SetBinContent( q+1, fitPim[p][q][x]->GetParameter(c));
						hFitPimQ[p][x][c]->SetBinError( q+1, fitPim[p][q][x]->GetParError(c));
					}

				}
			}

			for( int c = 0; c < 4; c++ ){
				//fitPipQ[p][x][c] = new TF1( Form("fitPipQ_%i_%i_%i",p,  x, c), "[0]*([1] - x)*([1] - x) + [2]");
				fitPipQ[p][x][c] = new TF1( Form("fitPipQ_%i_%i_%i",p,  x, c), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x");
				//fitPipQ[p][x][c]->SetParameters(10, 4, 1);
				fitPipQ[p][x][c]->SetParameters(10, 4, 1, 1);
				hFitPipQ[p][x][c]->Fit(Form("fitPipQ_%i_%i_%i",p, x, c));
				
				fitPipQ[p][x][c]->SetLineColor(kBlue);
				hFitPipQ[p][x][c]->SetLineColor(kAzure);;
				hFitPipQ[p][x][c]->Draw();
				fitPipQ[p][x][c]->Draw("SAME");

				//fitPimQ[p][x][c] = new TF1( Form("fitPimQ_%i_%i_%i",p,  x, c), "[0]*([1] - x)*([1] - x) + [2]");
				fitPimQ[p][x][c] = new TF1( Form("fitPimQ_%i_%i_%i",p,  x, c), "[0] + [1]*x + [2]*x*x + [3]*x*x*x");
				//fitPimQ[p][x][c]->SetParameters(10, 4, 1);
				fitPimQ[p][x][c]->SetParameters(10, 4, 1, 1);
				hFitPimQ[p][x][c]->Fit(Form("fitPimQ_%i_%i_%i",p, x, c));
				
				fitPimQ[p][x][c]->SetLineColor(kRed);
				hFitPimQ[p][x][c]->SetLineColor(kMagenta);;
				hFitPimQ[p][x][c]->Draw("SAME");
				fitPimQ[p][x][c]->Draw("SAME");




				canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
				canvas.Clear();
			}
		
			for( int i = 0; i < 4; i++ ){
				for( int j = 0; j < 3; j++ ){
					if( fitPipQ[p][x][i]->GetParError(j) != 0 && fitPipQ[p][x][i]->GetParameter(j) != 0){
						hFitPipX[p][i][j]->SetBinContent( x + 1, fitPipQ[p][x][i]->GetParameter(j) );	
						hFitPipX[p][i][j]->SetBinError( x + 1, fitPipQ[p][x][i]->GetParError(j) );	
					}
					//if( fitPipX_f[x][i]->GetParError(j) == 0 && fitPipX_f[x][i]->GetParameter(j) != 0){
					//	fitPipQ[i][j]->SetBinContent( x + 1, fitPipX_f[x][i]->GetParameter(j) );	
					//	fitPipQ[i][j]->SetBinError( x + 1, fitPipX_f[x][i]->GetParameter(j)*.1 );	
					//}
					//fitPimQ[i][j]->SetBinContent( x + 1, fitPimX_f[x][i]->GetParameter(j) );	
					if( fitPimQ[p][x][i]->GetParError(j) != 0 && fitPimQ[p][x][i]->GetParameter(j) != 0){
						hFitPimX[p][i][j]->SetBinContent( x + 1, fitPimQ[p][x][i]->GetParameter(j) );	
						hFitPimX[p][i][j]->SetBinError( x + 1, fitPimQ[p][x][i]->GetParError(j) );	
					}
				}
			}
		}
		for( int i = 0; i < 4; i++ ){
			for( int j = 0; j < 3; j++ ){
				fitPipX[p][i][j] = new TF1( Form("fitPipX_%i_%i_%i",p, i, j), "[0] + [1]*x + [2]*x*x + [3]*x*x*x");
				hFitPipX[p][i][j]->Fit( Form("fitPipX_%i_%i_%i", p, i, j) );
				fitPipX[p][i][j]->Write(); 
				
				fitPimX[p][i][j] = new TF1( Form("fitPimX_%i_%i_%i",p, i, j), "[0] + [1]*x + [2]*x*x + [3]*x*x*x");
				hFitPimX[p][i][j]->Fit( Form("fitPimX_%i_%i_%i", p, i, j) );
				fitPimX[p][i][j]->Write(); 
			
			
				//hFitPipX[p][i][j]->SetLineColor(kRed);
				//hFitPipX[p][i][j]->Draw("SAME");
				hFitPimX[p][i][j]->SetLineColor(kAzure);
				hFitPimX[p][i][j]->Draw("Same");
				canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
				canvas.Clear();
			}
		}
	}
	


	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
	outFile->Close();
}

