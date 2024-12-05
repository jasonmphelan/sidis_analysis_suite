//Adds new branch for final acceptance matching in 2d

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
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "cut_values.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "reader.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;
using namespace cutVals; 


int main( int argc, char** argv){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input Path] [Output File] [# of input files] [File Type] [Beam Energy]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\nnFiles " << atoi(argv[3]) << "\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
       	int runType = atoi(argv[3]);
       	double EBeam = atof(argv[4]);

	TFile * inFile = new TFile(in_name);
	
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);
	anal.loadMatchingFunctions();

	
	cout<<"Completed good event list... \n";
        //cout<<"Number of event (no acc) : "<< isGoodPion
	//Define out tree and files
	TTree * old_tree_upd = (TTree *) inFile->Get("ePi");
        TFile * outFile = new TFile(out_name, "RECREATE");

	cout<<"Created new tree and outfile\n";

	std::vector<pion> pi;
 	std::vector<bool> isGoodPion_fit;
	int EBeam_run = EBeam;
	
	old_tree_upd->SetBranchAddress("pi", &pi );

	TTree * outTree = old_tree_upd->CloneTree(0);
        TBranch * branch_1 = outTree->Branch("EBeam_run", &EBeam_run );
        TBranch * branch_2 = outTree->Branch("isGoodPion_fit", &isGoodPion_fit );
        
        

        int nEvents = old_tree_upd->GetEntries();

	cout<<"Starting event loop\n";

        for(int ev = 0; ev < nEvents; ev++){
		isGoodPion_fit.clear();
		pi.clear();
                if( ev %10000 == 0 ){cout <<ev<<" / "<<nEvents<<std::endl;}
				
                old_tree_upd->GetEntry(ev);

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			isGoodPion_fit.push_back(false);
                	if(!anal.applyPionKinematicCuts(pi[i])){ continue; }
			if ( anal.applyAcceptanceMatching(pi[i], 2) ){
				isGoodPion_fit[i] = true;
			}
		}

                outTree->Fill();
        }

        cout<<"Writing to file\n";
        outFile->cd();
        outTree->Write();
        outFile->Close();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout<<"Done. Elapsed Time : "<<elapsed.count()<<std::endl;

	return 0;
}
