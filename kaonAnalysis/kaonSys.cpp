#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TClass.h"
#include "TKey.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST


using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;


int main( int argc, char** argv){

	if( argc < 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [OutFile Base] [Infile bases...]\n";
		return -1;
	}

	double E_list[3] = {10.2, 10.4, 10.6};

	TString out_base = argv[1];
	std::vector<TString> in_base;

	if( argc == 2 ){
		in_base.push_back( argv[1] );
	}	
	else if( argc == 3 ){
		in_base.push_back( argv[2] );
	}
	else{
		for( int i = 2; i < argc; i++ ){
			in_base.push_back( argv[i] );
		}
	}	


	TFile * outFile = new TFile( out_base + "_ratio.root", "RECREATE" );
		
	TFile * inFile_1 = new TFile( out_base + ".root" );
	TFile * inFile_2 = new TFile( out_base + "_no_match.root" );

	TIter keyList(inFile_1->GetListOfKeys());
	TKey *key;
	outFile->cd();

  	while ((key = (TKey*)keyList())) {
      		TClass *cl = gROOT->GetClass(key->GetClassName());
      		if (cl->InheritsFrom("TH3")){
			TString keyName = (TString)key->GetName();
			cout<<"Doing : "<<keyName<<std::endl;
			TH3F * h1 = (TH3F*)inFile_1->Get(keyName);
			TH3F * h2 = (TH3F*)inFile_2->Get(keyName);


			h1->Divide(h2);

			h1->Write();
		}
    
	}
	outFile->Close();
}
