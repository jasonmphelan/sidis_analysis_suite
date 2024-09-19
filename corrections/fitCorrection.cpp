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

	if( argc <4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File (no extension)] [correction type ( 0=bin, 1=acc)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
	int corrType = atoi(argv[3]);

       	TFile * outFile = new TFile((TString) CORR_PATH + "/correctionFiles/" + out_name + ".root", "RECREATE");
       	TFile * inFile = new TFile((TString) CORR_PATH + "/correctionFiles/" + in_name + ".root");
	
	TString histName;
	if( corrType  == 0 ){
		histName = "hBinMigration";
	}
	else if ( corrType == 1 ){
		histName = "hAccCorrection";
	}
	else{
		cerr << "Incorrect correction type \n";
		return -1;
	}

	cout<<"GETTING HISTS\n";

	TH3F * kaonCorr_p = (TH3F *)inFile->Get( histName + "P" );
	TH3F * kaonCorr_m = (TH3F *)inFile->Get( histName + "M" );
	
	//Project to 1D histograms
	TH1F * kaonCorr_p_1d[bins_Q2][bins_xB];
	TH1F * kaonCorr_m_1d[bins_Q2][bins_xB];
	TF1 *fitPim[bins_Q2][bins_xB];
	TF1 *fitPip[bins_Q2][bins_xB];


	TH1F * fitPipX[bins_Q2][4];
	TH1F * fitPimX[bins_Q2][4];
	TF1 *fitPimX_f[bins_Q2][4];
	TF1 *fitPipX_f[bins_Q2][4];

	TH1F * fitPipQ[4][3];
	TH1F * fitPimQ[4][3];

	TString co_names[4] = {"a", "b", "c", "d"};

	for( int q = 0; q < bins_Q2; q++ ){
		for( int c = 0; c < 4; c++ ){
			fitPipX[q][c] = new TH1F(Form("fitPipX_%i_%i", q, c ), "", bins_xB, xB_min, xB_max);
			fitPimX[q][c] = new TH1F(Form("fitPimX_%i_%i", q, c ), "", bins_xB, xB_min, xB_max);

		}
	}
	for( int c = 0; c < 4; c++ ){
		for( int k = 0; k < 3; k++ ){
			fitPipQ[c][k] = new TH1F( Form("fitPipQ_%i_%i", c, k), "", bins_Q2, Q2_min, Q2_max );
			fitPimQ[c][k] = new TH1F( Form("fitPimQ_%i_%i", c, k), "", bins_Q2, Q2_min, Q2_max );
		}
	}

	outFile->cd();

	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int q = 0; q <  bins_Q2; q++ ){
		for( int x = 0; x < bins_xB; x++ ){
			kaonCorr_p_1d[q][x] = new TH1F( histName + Form("_p_1d_%i_%i", q, x), ";z;w+", bins_Z, Z_min, Z_max );
			kaonCorr_m_1d[q][x] = new TH1F( histName + Form("_m_1d_%i_%i", q, x), ";z;w-", bins_Z, Z_min, Z_max );
			
			

			double p_fit_min = 0;
			double p_fit_max = 1;
			double m_fit_min = 0;
			double m_fit_max = 1;

			for( int z = 0; z < bins_Z; z++ ){
				if( !isnan(kaonCorr_p->GetBinError( x+1, q+1, z+1) )  ){
				//if( !isnan(kaonCorr_p->GetBinContent( x+1, q+1, z+1) ) || isfinite(kaonCorr_p->GetBinContent( x+1, q+1, z+1) )  ){
					kaonCorr_p_1d[q][x]->SetBinContent( z+1, kaonCorr_p->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_p_1d[q][x]->SetBinError( z+1, kaonCorr_p->GetBinError( x+1, q+1, z+1 ));
				
				}
				else{
					kaonCorr_p_1d[q][x]->SetBinContent( z+1, 0);
				}
					
				//if( !isnan(kaonCorr_m->GetBinContent( x+1, q+1, z+1) ) || isfinite(kaonCorr_m->GetBinContent( x+1, q+1, z+1) )  ){
				if( !isnan(kaonCorr_m->GetBinError( x+1, q+1, z+1) )  ){
					kaonCorr_m_1d[q][x]->SetBinContent( z+1, kaonCorr_m->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_m_1d[q][x]->SetBinError( z+1, kaonCorr_m->GetBinError( x+1, q+1, z+1 ));
				
				}
				else{
					kaonCorr_m_1d[q][x]->SetBinContent( z+1, 0);
				}
			}
			for( int z = 0; z < bins_Z; z++ ){
				if( p_fit_min == 0 &&  kaonCorr_p->GetBinContent( x+1, q+1, z+1 ) > 0 ){
					p_fit_min = kaonCorr_p_1d[q][x]->GetBinCenter(z+1) - .025;
				}
				if( p_fit_max == 1 && p_fit_min > 0 &&  kaonCorr_p->GetBinContent( x+1, q+1, z+1 ) == 0){
					p_fit_max = kaonCorr_p_1d[q][x]->GetBinCenter(z+1) + .025;
				}
				if( m_fit_min == 0 &&  kaonCorr_m->GetBinContent( x+1, q+1, z+1 ) > 0 ){
					m_fit_min = kaonCorr_m_1d[q][x]->GetBinCenter(z+1) - .025;
				}
				if( m_fit_max == 1 && m_fit_min > 0 &&  kaonCorr_m->GetBinContent( x+1, q+1, z+1 ) == 0){
					m_fit_max = kaonCorr_m_1d[q][x]->GetBinCenter(z+1) + .025;
				}
			}

			fitPip[q][x] = new TF1( Form("fitPip_%i_%i", q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", p_fit_min, p_fit_max );
			//fitPip[q][x] = new TF1( Form("fitPip_%i_%i", q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
			fitPip[q][x]->SetParameters( 1, 1, 1, 1 );
			
			fitPim[q][x] = new TF1( Form("fitPim_%i_%i", q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", m_fit_min, m_fit_max );
			//fitPim[q][x] = new TF1( Form("fitPim_%i_%i", q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
			fitPim[q][x]->SetParameters( 1, 1, 1, 1 );
			
			kaonCorr_p_1d[q][x]->Fit( Form("fitPip_%i_%i", q, x) );
			kaonCorr_m_1d[q][x]->Fit( Form("fitPim_%i_%i", q, x) );
			kaonCorr_p_1d[q][x]->Draw();
			//fitPip[q][x]->Draw("SAME");
			//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
			//canvas.Clear();
			
			kaonCorr_p_1d[q][x]->Write();
			kaonCorr_m_1d[q][x]->Write();
		
		

			for( int c = 0; c < 4; c++ ){
				if( fitPip[q][x]->GetParError(c) != 0 ){
					fitPipX[q][c]->SetBinContent( x+1, fitPip[q][x]->GetParameter(c));
					fitPipX[q][c]->SetBinContent( x+1, fitPip[q][x]->GetParError(c));
				}
				if( fitPim[q][x]->GetParError(c) != 0 ){
					fitPimX[q][c]->SetBinContent( x+1, fitPim[q][x]->GetParameter(c));
					fitPimX[q][c]->SetBinContent( x+1, fitPim[q][x]->GetParError(c));
				}
			}
		}
		fitPipX_f[q][0] = new TF1( Form("fitPipX_f_%i_%i", q, 0), "[0]*([1] - x)*([1] - x) + [2]");
		fitPipX_f[q][0]->SetParameters(1, 1, 1);
		fitPipX[q][0]->SetLineColor(kRed);
		//fitPipX[q][0]->Fit("chebyshev2");
		fitPipX[q][0]->Fit(Form("fitPipX_f_%i_%i", q, 0));
		//fitPipX[q][0]->Draw("E1");
		//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
		//canvas.Clear();
		
		fitPipX_f[q][1] = new TF1( Form("fitPipX_f_%i_%i", q, 1), "[0]*([1] - x)*([1] - x) + [2]");
		fitPipX_f[q][1]->SetParameters(1, 1, 1);
		fitPipX[q][1]->SetLineColor(kAzure);
		//fitPipX[q][1]->Fit("chebyshev2");
		fitPipX[q][1]->Fit(Form("fitPipX_f_%i_%i", q, 1));
		//fitPipX[q][1]->Draw("E1");
		//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
		//canvas.Clear();
		
		fitPipX_f[q][2] = new TF1( Form("fitPipX_f_%i_%i", q, 2), "[0]*([1] - x)*([1] - x) + [2]");
		fitPipX_f[q][2]->SetParameters(1, 1, 1);
		fitPipX[q][2]->SetLineColor(kGreen);
		//fitPipX[q][2]->Fit("chebyshev2");
		fitPipX[q][2]->Fit(Form("fitPipX_f_%i_%i", q, 2));
		//fitPipX[q][2]->Draw("E1");
		//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
		//canvas.Clear();
		
		fitPipX_f[q][3] = new TF1( Form("fitPipX_f_%i_%i", q, 3), "[0]*([1] - x)*([1] - x) + [2]");
		fitPipX_f[q][3]->SetParameters(1, 1, 1);
		//fitPipX[q][3]->Fit("chebyshev2");
		//fitPipX[q][3]->SetLineColor(kBlack);
		fitPipX[q][3]->Fit(Form("fitPipX_f_%i_%i", q, 3));
		//fitPipX[q][3]->Draw("E1");
		//canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
		//canvas.Clear();

		for( int i = 0; i < 4; i++ ){
			for( int j = 0; j < 3; j++ ){
				cout<<" BIN VAL : "<<fitPipX_f[q][i]->GetParameter(j)<<std::endl;
				if( fitPipX_f[q][i]->GetNDF() != 0){
					fitPipQ[i][j]->SetBinContent( q + 1, fitPipX_f[q][i]->GetParameter(j) );	
				}
				//fitPimQ[i][j]->SetBinContent( q + 1, fitPimX_f[q][i]->GetParameter(j) );	
			}
		}

	}
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j < 3; j++ ){
			fitPipQ[i][j]->SetLineColor(kRed);
			fitPipQ[i][j]->Draw();
			//fitPimQ[i][j]->SetLineColor(kBlue);
			//fitPimQ[i][j]->Draw();
			canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
			canvas.Clear();
		}
	}
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");

	outFile->Close();
}

