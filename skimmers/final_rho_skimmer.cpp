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
#include "clashit.h"


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

const TLorentzVector p_rest( 0, 0, 0, 0.938 );

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File] [Cut Level (0, 1)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\nCut Level " << atoi(argv[3]) << "\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];


        TFile * file_rec = new TFile(in_name, "UPDATE");
	//SIDISatBAND_auxiliary aux(1); //For kinematical cuts
    	//aux.loadCutValues("macros/cuts/BANDcutValues.csv", -1);

	//Load Match Functions
	TFile * matchFile = new TFile("/work/clas12/users/jphelan/SIDIS_at_BAND/matchCut.root", "UPDATE");
	TF1 * max[6][2];
	TF1 * min[6][2];
	TString chargeType[2] = {"pip", "pim"};
	
	for( int i = 0; i < 6; i++ ){
		
		for( int j = 0; j < 2; j++ ){
			max[i][j] = (TF1 *)matchFile->Get(Form("max_%i_", i) + chargeType[j]);
			min[i][j] = (TF1 *)matchFile->Get(Form("min_%i_", i) + chargeType[j]);
		}
	}

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<double> Ebeam_ptr(reader_rec, "Ebeam");
        TTreeReaderValue<TLorentzVector> p_e(reader_rec, "p_e");
        TTreeReaderValue<TLorentzVector> q(reader_rec, "q");
        TTreeReaderValue<double> Q2_ptr(reader_rec, "q2");
        TTreeReaderValue<double> xB_ptr(reader_rec, "xB");
        TTreeReaderValue<double> W_ptr(reader_rec, "w");
        TTreeReaderValue<double> y_ptr(reader_rec, "y");
        TTreeReaderValue<double> omega_ptr(reader_rec, "omega");

        TTreeReaderArray<TLorentzVector> p_pi(reader_rec, "p_pi");
        TTreeReaderArray<clashit> pi(reader_rec, "pi");
        TTreeReaderArray<double> Z(reader_rec, "Z");
        TTreeReaderArray<double> M_x(reader_rec, "M_x");
        //TTreeReaderArray<int> sector(reader_rec, "pi_sector_DC");
	TTreeReaderArray<int> charge(reader_rec, "charge");

	//Define good event list and additional variables for output branches
        TEventList * good_events = new TEventList();
        
	std::vector<double> M_rho;
	std::vector<double> M_x_rho;
	std::vector<int> leadIdx_vec;
	std::vector<std::vector<bool>> isGoodPion_temp;

	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		//Skim for potential neutral rho events
		if( (int) (p_pi.end() - p_pi.begin()) != 2 ){continue;} 
		if( charge[0] == charge[1] ){ continue; }

		double Ebeam = *Ebeam_ptr;
                double Q2 = *Q2_ptr;
		double W = *W_ptr;
                double xB = * xB_ptr;
                double y = *y_ptr;
		double omega = * omega_ptr;
	
		//require good electron
		if( W < 2.5 ) { continue; }
		if( Q2 < Q2_min || Q2 > Q2_max ) { continue; }
                if( xB < xB_min || xB > xB_max ) { continue; }
                if( y > 0.75 ) { continue; }
                if( p_e->Vect().Mag() <3. || p_e->Vect().Mag() > 10.6 ) { continue; }
                if( p_e->Theta()*rad_to_deg < 5. || p_e->Theta()*rad_to_deg > 35. ){ continue; }

		double lead_Z_temp = -1;
		int lead_idx_temp = -1;
		
		std::vector<bool> isGoodPion_event;
		for( int i = 0; i < (int) ( p_pi.end() - p_pi.begin() ); i++ ){
			isGoodPion_event.push_back(false);
               
                        int sector_i = pi[i].getDC_sector();
			
			//Pion cuts
                        if ( ( M_x[i] < 1.7 || M_x[i] > 5.) ) { continue; }
                        if ( p_pi[i].Vect().Mag() < 1.25 || p_pi[i].Vect().Mag() > 5. ) { continue; }
                        if ( Z[i] < Z_min  ||  Z[i] > Z_max ) { continue; }
                        if ( p_pi[i].Theta()*rad_to_deg < 5. || p_pi[i].Theta()*rad_to_deg > 35 ){ continue; }
                	
			double theta = p_pi[i].Theta()*rad_to_deg;

			double acc_map_pip_min = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/p_pi[i].Vect().Mag();                      
                        double acc_map_pim_min = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/p_pi[i].Vect().Mag();
		
			//acceptance matching
			if ( theta > acc_map_pip_min && theta > acc_map_pim_min ){
				if( Z[i] > lead_Z_temp ){
					lead_Z_temp = Z[i];
					lead_idx_temp = i;
				}	
				//if pass all cuts, save as good pion
				isGoodPion_event[i] = true;
			}

		}
		
		if(lead_idx_temp > -1){
			good_events->Enter(event_count);
			
			isGoodPion_temp.push_back(isGoodPion_event);
			leadIdx_vec.push_back( lead_idx_temp );
			M_rho.push_back( (p_pi[0] + p_pi[1]).Mag() );
			M_x_rho.push_back( (p_rest + *q - p_pi[0] - p_pi[1]).Mag() );
		}
	}
	
	cout<<"Completed good event list... \n";
	//Define out tree and files
	TTree * old_tree_upd = (TTree *) reader_rec.GetTree();
        TFile * outFile = new TFile(out_name, "RECREATE");

	cout<<"Created new tree and outfile\n";

        int nLeadIdx;
	double M_x_rho_event, M_rho_event;

 	std::vector<bool> isGoodPion;

	cout<<"Declared new variables\n";

	TTree * outTree = old_tree_upd->CloneTree(0);

	outTree->Print();        

	TBranch * branch = outTree->Branch("nLeadIdx", &nLeadIdx, "nLeadIdx/I" );
        TBranch * branch_2 = outTree->Branch("isGoodPion", &isGoodPion );
        TBranch * branch_3 = outTree->Branch("M_x_rho", &M_x_rho_event );
        TBranch * branch_4 = outTree->Branch("M_rho", &M_rho_event );
        
        outTree->Print();

        int nEvents = good_events->GetN();

	cout<<"Starting event loop\n";

        for(int ev = 0; ev < nEvents; ev++){

                if( ev %10000 == 0 ){cout <<ev<<" / "<<nEvents<<std::endl;}
				
                int entry  = good_events->GetEntry(ev);
                old_tree_upd->GetEntry(entry);
                
		nLeadIdx = leadIdx_vec[ev];
		M_x_rho_event = M_x_rho[ev];
		M_rho_event = M_rho[ev];

		isGoodPion = isGoodPion_temp[ev];

                outTree->Fill();
        }

        cout<<"Writing to file\n";
        outFile->cd();
	//accCharge->Write("accCharge");
        outTree->Write();
        outFile->Close();
	return 0;
}
