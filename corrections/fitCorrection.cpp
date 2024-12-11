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
	else if ( corrType == 2 ){
		histName = "hMcCorrection";
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


	outFile->cd();

	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int x = 0; x < bins_xB; x++ ){
		for( int q = 0; q <  bins_Q2; q++ ){
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
				//else{
				//	kaonCorr_p_1d[q][x]->SetBinContent( z+1, 0);
				//}
					
				//if( !isnan(kaonCorr_m->GetBinContent( x+1, q+1, z+1) ) || isfinite(kaonCorr_m->GetBinContent( x+1, q+1, z+1) )  ){
				if( !isnan(kaonCorr_m->GetBinError( x+1, q+1, z+1) )   ){
					kaonCorr_m_1d[q][x]->SetBinContent( z+1, kaonCorr_m->GetBinContent( x+1, q+1, z+1 ));
					kaonCorr_m_1d[q][x]->SetBinError( z+1, kaonCorr_m->GetBinError( x+1, q+1, z+1 ));
				
				}
				//else{
				//	kaonCorr_m_1d[q][x]->SetBinContent( z+1, 0);
				//}
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
			
			
			cout<<"##############################################################\n\n"<<kaonCorr_p_1d[q][x]->GetEntries() << std::endl;
			fitPip[q][x] = new TF1( Form("fitPip_%i_%i", q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", p_fit_min, p_fit_max );
			//fitPip[q][x] = new TF1( Form("fitPip_%i_%i", q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
			fitPip[q][x]->SetParameters( 1, 1, 1, 1 );
			
			fitPim[q][x] = new TF1( Form("fitPim_%i_%i", q, x), "[0]+ [1]*x + [2]*x*x + [3]*x*x*x", m_fit_min, m_fit_max );
			//fitPim[q][x] = new TF1( Form("fitPim_%i_%i", q, x), "[0]+ [1]*pow((1. - x), [2])", Z_min, Z_max );
			fitPim[q][x]->SetParameters( 1, 1, 1, 1 );
			
			if( kaonCorr_p_1d[q][x]->GetEntries() == 0 ){
				fitPip[q][x]->SetParameters( 0, 0, 0, 0 );
			}
			else{
				kaonCorr_p_1d[q][x]->Fit( Form("fitPip_%i_%i", q, x) );
			}
			if( kaonCorr_m_1d[q][x]->GetEntries() == 0 ){
				fitPim[q][x]->SetParameters( 0, 0, 0, 0);
			}
			else{
				kaonCorr_m_1d[q][x]->Fit( Form("fitPim_%i_%i", q, x) );
			}
			kaonCorr_p_1d[q][x]->Draw();
			fitPip[q][x]->Draw("SAME E3");
			fitPip[q][x]->Write();
			canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
			canvas.Clear();
			kaonCorr_m_1d[q][x]->Draw();
			fitPim[q][x]->Draw("same E3");
			canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
			fitPim[q][x]->Write();
			canvas.Clear();
			
		}
		canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");	
		canvas.Clear();
		cout<<"FINISHED FITS\n";
		

	}
	cout<<"Finished Fits\n";
	
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");

	outFile->Close();
}

