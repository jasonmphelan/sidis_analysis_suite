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
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
#include "reader.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#endif

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;
using namespace cutVals; 


int main( int argc, char** argv){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 6 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input Path] [Output File] [# of input files] [File Type] [Beam Energy]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\nnFiles " << atoi(argv[3]) << "\n";
	TString in_name = argv[1];
       	TString out_name = argv[2];
       	int nFiles = atoi(argv[3]);
       	int runType = atoi(argv[4]);
       	double EBeam = atof(argv[5]);

	reader skimReader;
	skimReader.setNumFiles( nFiles);
	skimReader.setRunType( runType );
	skimReader.setEnergy( EBeam );

	TChain * chain = new TChain("ePi");
	skimReader.getRunSkimsByName(chain, in_name);

        //TFile * file_rec = new TFile(in_name, "UPDATE");
	
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions();
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));


	//Load input tree
        
        //TTreeReader reader_rec( chain );
	//TTreeReaderValue<electron> e(reader_rec, "e");
        //TTreeReaderArray<pion> pi(reader_rec, "pi");
		
		// Pointers to branch objects
	electron* e = nullptr;
	std::vector<pion>* pi = nullptr;

	genElectron* e_gen = nullptr;
	std::vector<genPion>* pi_gen = nullptr;

	chain->SetBranchAddress("e", &e);
	chain->SetBranchAddress("pi", &pi);
	if( runType == 1 ){
		chain->SetBranchAddress("e_gen", &e_gen);
		chain->SetBranchAddress("pi_gen", &pi_gen);
	}
	double event_total = chain->GetEntries();

	//Define good event list and additional variables for output branches
    TEventList * good_events = new TEventList();
        
	//std::vector<int> leadIdx_vec;
	//std::vector<int> leadIdx_no_acc_vec;
	//std::vector<int> leadIdx_3d_vec;
	std::vector<std::vector<bool>> isGoodPion_temp;
	std::vector<std::vector<bool>> isGoodGenPion_temp;
	std::vector<std::vector<bool>> isGoodPion_no_acc_temp;
	std::vector<std::vector<bool>> isGoodPion_3d_temp;

	//int event_total = reader_rec.GetEntries();

	for ( int event_count = 0; event_count < event_total; event_count++ ){
		chain->GetEntry(event_count);
	//while (reader_rec.Next()) {
    	//int event_count = reader_rec.GetCurrentEntry();

		if(event_count%10000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		
		std::vector<bool> isGoodPion_event;
		std::vector<bool> isGoodPion_no_acc_event;
		std::vector<bool> isGoodPion_3d_event;
		std::vector<bool> isGoodGenPion_event;
		bool e_pass = anal.applyElectronKinematicCuts( *e );
		bool e_gen_pass = true;
		if( runType == 1 ) e_gen_pass = anal.applyElectronKinematicCuts( *e_gen );

		if( !e_pass && !e_gen_pass ){ continue; }

		for( int i = 0; i < (int) ( pi->end() - pi->begin() ); i++ ){
			isGoodPion_event.push_back(false);
			isGoodGenPion_event.push_back(false);
			isGoodPion_no_acc_event.push_back(false);
			isGoodPion_3d_event.push_back(false);
		     
			bool pi_pass = anal.applyPionKinematicCuts((*pi)[i]);
			bool pi_gen_pass = true;
			if(runType == 1 ) pi_gen_pass = anal.applyPionKinematicCuts((*pi_gen)[i]);

			if( !pi_pass && !pi_gen_pass ) continue;

			isGoodPion_no_acc_event[i] = (pi_pass && e_pass);
			isGoodGenPion_event[i] = (pi_gen_pass && e_gen_pass);
			
			double p_pi = (*pi)[i].get3Momentum().Mag();
			double theta_pi = (*pi)[i].get3Momentum().Theta();
			double phi_pi = (*pi)[i].get3Momentum().Phi();

			if(anal.applyAcceptanceMap( e->get3Momentum().Mag(), rad_to_deg*e->get3Momentum().Phi(), rad_to_deg*e->get3Momentum().Theta(), 1 ) <0) continue;
			if(anal.applyAcceptanceMap( p_pi, rad_to_deg*(*pi)[i].get3Momentum().Phi(), rad_to_deg*(*pi)[i].get3Momentum().Theta(), 1 ) < 0 ) continue;
			
			if ( anal.applyAcceptanceMatching((*pi)[i], 2) ){
				isGoodPion_event[i] = true;
			}
			

			if ( anal.applyAcceptanceMap( p_pi, rad_to_deg*(*pi)[i].get3Momentum().Phi(), rad_to_deg*(*pi)[i].get3Momentum().Theta(), 1 ) >= 0 &&
							anal.applyAcceptanceMap( p_pi, rad_to_deg*(*pi)[i].get3Momentum().Phi(), rad_to_deg*(*pi)[i].get3Momentum().Theta(), 2 ) >= 0 ){
				
		
				isGoodPion_3d_event[i] = true;
			}
			

		}

		bool hasGoodReco = std::any_of(isGoodPion_no_acc_event.begin(),
                                   isGoodPion_no_acc_event.end(),
                                   [](bool v){ return v; });

    	bool hasGoodGen  = std::any_of(isGoodGenPion_event.begin(),
                                   isGoodGenPion_event.end(),
                                   [](bool v){ return v; });
		
		if(hasGoodReco || hasGoodGen){
			good_events->Enter(event_count);
			isGoodGenPion_temp.push_back(isGoodGenPion_event);	
			isGoodPion_no_acc_temp.push_back(isGoodPion_no_acc_event);
			isGoodPion_temp.push_back(isGoodPion_event);
			isGoodPion_3d_temp.push_back(isGoodPion_3d_event);
		}
	}
	
	cout<<"Completed good event list... \n";
    //cout<<"Number of event (no acc) : "<< isGoodPion
	//Define out tree and files
	TTree * old_tree_upd = (TTree *) chain->GetTree();
    TFile * outFile = new TFile(out_name, "RECREATE");

	cout<<"Created new tree and outfile\n";


 	std::vector<bool> isGoodPion;
 	std::vector<bool> isGoodPion_no_acc;
 	std::vector<bool> isGoodPion_3d;
 	std::vector<bool> isGoodGenPion;


	TTree * outTree = old_tree_upd->CloneTree(0);

    TBranch * branch_1 = outTree->Branch("isGoodPion", &isGoodPion );
    TBranch * branch_2 = outTree->Branch("isGoodGenPion", &isGoodGenPion );
        
    TBranch * branch_3 = outTree->Branch("isGoodPion_no_acc", &isGoodPion_no_acc);
        
    TBranch * branch_4 = outTree->Branch("isGoodPion_3d", &isGoodPion_3d);

    int nEvents = good_events->GetN();

	cout<<"Starting event loop\n";

        for(int ev = 0; ev < nEvents; ev++){

            if( ev %10000 == 0 ){cout <<ev<<" / "<<nEvents<<std::endl;}
				
            int entry  = good_events->GetEntry(ev);
            old_tree_upd->GetEntry(entry);
		
		isGoodGenPion = isGoodGenPion_temp[ev];
		isGoodPion_no_acc = isGoodPion_no_acc_temp[ev]; 
		isGoodPion_3d = isGoodPion_3d_temp[ev]; 
		isGoodPion = isGoodPion_temp[ev];
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
