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


#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<clashit>+;
#pragma link C++ class vector<clashit>+;
#endif

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;
using namespace cutVals; 


int main( int argc, char** argv){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc <4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File] [Cut Level (0, 1)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\nCut Level " << atoi(argv[3]) << "\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
       	int cut_type = atoi(argv[3]);


        TFile * file_rec = new TFile(in_name, "UPDATE");
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");

        TTreeReaderArray<pion> pi(reader_rec, "pi");

	//Define good event list and additional variables for output branches
        TEventList * good_events = new TEventList();
        
	std::vector<int> leadIdx_vec;
	std::vector<int> leadIdx_no_acc_vec;
	//std::vector<int> leadIdx_3d_vec;
	std::vector<std::vector<bool>> isGoodPion_temp;
	std::vector<std::vector<bool>> isGoodPion_no_acc_temp;
	//std::vector<std::vector<bool>> isGoodPion_3d_temp;

	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}


		double lead_Z_temp = -1;
		double lead_Z_no_acc_temp = -1;               
		//double lead_Z_3d_temp = -1;               
 
		int lead_idx_temp = -1;
                int lead_idx_no_acc_temp = -1;
                //int lead_idx_3d_temp = -1;
		
		std::vector<bool> isGoodPion_event;
		std::vector<bool> isGoodPion_no_acc_event;
		//std::vector<bool> isGoodPion_3d_event;

		if( !anal.applyElectronKinematicCuts( *e ) ){ continue; }

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			isGoodPion_event.push_back(false);
			isGoodPion_no_acc_event.push_back(false);;
			//isGoodPion_3d_event.push_back(false);;
               
                	if(!anal.applyPionKinematicCuts(pi[i])){ continue; }

			if( pi[i].getZ() > lead_Z_no_acc_temp ){
				lead_Z_no_acc_temp = pi[i].getZ();
				lead_idx_no_acc_temp = i;
			}	
        		

			isGoodPion_no_acc_event[i] = true;;

			if ( anal.applyAcceptanceMatching(pi[i], 2) ){
				if( pi[i].getZ() > lead_Z_temp ){
					lead_Z_temp = pi[i].getZ();
					lead_idx_temp = i;
				}	
		
				isGoodPion_event[i] = true;
			}
			/*
			if( acceptance_match_3d_cont( phi, theta, p, match3d ) ){
				if( Z[i] > lead_Z_3d_temp ){
					lead_Z_3d_temp = Z[i];
					lead_idx_3d_temp = i;
				}	
		
				isGoodPion_3d_event[i] = true;
			}
			*/

		}
		
		if(lead_idx_no_acc_temp > -1){
			good_events->Enter(event_count);
	
			isGoodPion_no_acc_temp.push_back(isGoodPion_no_acc_event);
			isGoodPion_temp.push_back(isGoodPion_event);
			//isGoodPion_3d_temp.push_back(isGoodPion_3d_event);

			leadIdx_no_acc_vec.push_back( lead_idx_no_acc_temp );
			leadIdx_vec.push_back( lead_idx_temp );
			//leadIdx_3d_vec.push_back( lead_idx_3d_temp );
		}
	}
	
	cout<<"Completed good event list... \n";
        //cout<<"Number of event (no acc) : "<< isGoodPion
	//Define out tree and files
	TTree * old_tree_upd = (TTree *) reader_rec.GetTree();
        TFile * outFile = new TFile(out_name, "RECREATE");

	cout<<"Created new tree and outfile\n";

        int nLeadIdx;
        int nLeadIdx_no_acc;
        //int nLeadIdx_3d;

 	std::vector<bool> isGoodPion;
 	std::vector<bool> isGoodPion_no_acc;
 	//std::vector<bool> isGoodPion_3d;


	TTree * outTree = old_tree_upd->CloneTree(0);

	TBranch * branch = outTree->Branch("nLeadIdx", &nLeadIdx, "nLeadIdx/I" );
        TBranch * branch_2 = outTree->Branch("isGoodPion", &isGoodPion );
        
	TBranch * branch_3 = outTree->Branch("nLeadIdx_no_acc", &nLeadIdx_no_acc, "nLeadIdx_no_acc/I" );
        TBranch * branch_4 = outTree->Branch("isGoodPion_no_acc", &isGoodPion_no_acc);
        
	//TBranch * branch_5 = outTree->Branch("nLeadIdx_3d", &nLeadIdx_3d, "nLeadIdx_3d/I" );
        //TBranch * branch_6 = outTree->Branch("isGoodPion_3d", &isGoodPion_3d);

        int nEvents = good_events->GetN();

	cout<<"Starting event loop\n";

        for(int ev = 0; ev < nEvents; ev++){

                if( ev %10000 == 0 ){cout <<ev<<" / "<<nEvents<<std::endl;}
				
                int entry  = good_events->GetEntry(ev);
                old_tree_upd->GetEntry(entry);
                
		nLeadIdx = leadIdx_vec[ev];
                nLeadIdx_no_acc = leadIdx_no_acc_vec[ev];
                //nLeadIdx_3d = leadIdx_3d_vec[ev];

		isGoodPion_no_acc = isGoodPion_no_acc_temp[ev]; 
		//isGoodPion_3d = isGoodPion_3d_temp[ev]; 
		isGoodPion = isGoodPion_temp[ev];

                outTree->Fill();
        }

        cout<<"Writing to file\n";
        outFile->cd();
	//accCharge->Write("accCharge");
        outTree->Write();
        outFile->Close();
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout<<"Done. Elapsed Time : "<<elapsed.count()<<std::endl;

	return 0;
}
