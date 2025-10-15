#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TObjArray.h"
#include "TObjString.h"
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

	if( argc < 1 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [OutFile Base] [Infile bases...]\n";
		return -1;
	}

	TString base = argv[1];

	TChain * ch = new TChain("ePi");
	ch->Add(base+"*.root");

	TObjArray* files = ch->GetListOfFiles();
	int nFiles = files->GetSize();

	if( nFiles != 100 ){
		cout<<"Missing Files... only found "<<nFiles<<endl;

		files->Print();
		return -1;
	}

	TFile * outFile = new TFile( base + ".root", "RECREATE" );
	TTree * outTree = ch->CloneTree(0);

	int nevent = ch->GetEntries();

	for ( int i = 0; i < nevent; i++ ){
		ch->GetEntry(i);
		//outTree->Fill();
	}

	outFile->cd();
	outTree->Write();
	outFile->Close();
}

