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
using std::isnan;
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

	TH3F * hCorr[2];
	hCorr[0] = (TH3F *)inFile->Get( histName + "P" );
	hCorr[1] = (TH3F *)inFile->Get( histName + "M" );
	
	TH3F * hSmooth[2];
	for( int i = 0; i < 2; i++ ){
		hSmooth[i] = (TH3F *)hCorr[i]->Clone();
	}
	

	//Project to 1D histograms

	outFile->cd();
	for( int ch = 0; ch < 2; ch++){
		for( int x = 0; x < bins_xB; x++ ){
			for( int q = 0; q <  bins_Q2; q++ ){
			
				TH1F * hCorr_1d = new TH1F( histName + Form("_p_1d_%i_%i", q, x), ";z;w+", bins_Z, Z_min, Z_max );

				//Fill 1d histogram

				for( int z = 0; z < bins_Z; z++ ){
					if( !isnan(hCorr[ch]->GetBinError( x+1, q+1, z+1) )  ){
						hCorr_1d->SetBinContent( z+1, hCorr[ch]->GetBinContent( x+1, q+1, z+1 ));
						hCorr_1d->SetBinError( z+1, hCorr[ch]->GetBinError( x+1, q+1, z+1 ));
					}
					else{
						hCorr_1d->SetBinContent( z+1, 0);
						hCorr_1d->SetBinError( z+1, 0);
					}
				}

				double min = hCorr_1d->GetBinLowEdge( hCorr_1d->FindFirstBinAbove(0) );
				double max = hCorr_1d->GetBinLowEdge( hCorr_1d->FindLastBinAbove(0) + 1 );
				
				if( hCorr_1d->Integral()>0){
					hCorr_1d->GetXaxis()->SetRangeUser( min, max);
					hCorr_1d->Smooth(1, "R");
				}

				for( int z = 0; z < bins_Z; z++){
					hSmooth[ch]->SetBinContent( x+1, q+1, z+1, hCorr_1d->GetBinContent( z + 1 ) );
				}
			}
			
		}
	
		hSmooth[ch]->Write();
	}
	

	outFile->Close();
}


