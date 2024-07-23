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

bool acceptance_match_3d( double phi_part, double theta, double p, int charge );
bool acceptance_match_3d_cont( double phi_part, double theta, double p, TF1 * fitFuncs[6][3]);

int main( int argc, char** argv){

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
	//SIDISatBAND_auxiliary aux(1); //For kinematical cuts
    	//aux.loadCutValues("macros/cuts/BANDcutValues.csv", -1);

	//Load Match Functions
	TFile * matchFile = new TFile("/work/clas12/users/jphelan/SIDIS_at_BAND/matchCut.root", "UPDATE");
	TF1 * max[6][2];
	TF1 * min[6][2];
	TString chargeType[2] = {"pip", "pim"};
	
	TFile * matchFile3D = new TFile("/work/clas12/users/jphelan/SIDIS_at_BAND/matchCut3D.root", "UPDATE");
	TF1 * match3d[6][3];

	for( int i = 0; i < 6; i++ ){
		
		match3d[i][0] = (TF1 *)matchFile3D->Get(Form("fTheta0_%i", i));
		match3d[i][1] = (TF1 *)matchFile3D->Get(Form("fPhi0_pip_%i", i));
		match3d[i][2] = (TF1 *)matchFile3D->Get(Form("fPhi0_pim_%i", i));

		for( int j = 0; j < 2; j++ ){
			max[i][j] = (TF1 *)matchFile->Get(Form("max_%i_", i) + chargeType[j]);
			min[i][j] = (TF1 *)matchFile->Get(Form("min_%i_", i) + chargeType[j]);
		}
	}

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<double> Ebeam_ptr(reader_rec, "Ebeam");
        TTreeReaderValue<TLorentzVector> p_e(reader_rec, "p_e");
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

        //TTreeReaderArray<bool> isMatched(reader_rec, "isMatched");
	
	//Define good event list and additional variables for output branches
        TEventList * good_events = new TEventList();
        
	std::vector<int> leadIdx_vec;
	std::vector<int> leadIdx_no_acc_vec;
	std::vector<int> leadIdx_3d_vec;
	std::vector<std::vector<bool>> isGoodPion_temp;
	std::vector<std::vector<bool>> isGoodPion_no_acc_temp;
	std::vector<std::vector<bool>> isGoodPion_3d_temp;

	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

		double Ebeam = *Ebeam_ptr;
                double Q2 = *Q2_ptr;
		double W = *W_ptr;
                double xB = * xB_ptr;
                double y = *y_ptr;
		double omega = * omega_ptr;
	
		if( W < 2.5 ) { continue; }
		if( Q2 < Q2_min || Q2 > Q2_max ) { continue; }
                if( xB < xB_min || xB > xB_max ) { continue; }
                if( y > 0.75 ) { continue; }
                if( p_e->Vect().Mag() <3. || p_e->Vect().Mag() > 10.6 ) { continue; }
                if( p_e->Theta()*rad_to_deg < 5. || p_e->Theta()*rad_to_deg > 35. ){ continue; }

		double lead_Z_temp = -1;
		double lead_Z_no_acc_temp = -1;               
		double lead_Z_3d_temp = -1;               
 
		int lead_idx_temp = -1;
                int lead_idx_no_acc_temp = -1;
                int lead_idx_3d_temp = -1;
		
		std::vector<bool> isGoodPion_event;
		std::vector<bool> isGoodPion_no_acc_event;
		std::vector<bool> isGoodPion_3d_event;

		for( int i = 0; i < (int) ( p_pi.end() - p_pi.begin() ); i++ ){
			isGoodPion_event.push_back(false);
			isGoodPion_no_acc_event.push_back(false);;
			isGoodPion_3d_event.push_back(false);;
               
			//cout<<"Found pion...\n";
         
                        int sector_i = pi[i].getDC_sector();
    			//if ( sector_i == 0 ) { continue;}			
			
                        if ( cut_type == 0 && ( M_x[i] < 1.7 || M_x[i] > 5.) ) { continue; }
                        if ( cut_type == 1 && ( M_x[i] < 1.5 || M_x[i] > 5.) ) { continue; }
                        if ( p_pi[i].Vect().Mag() < 1.25 || p_pi[i].Vect().Mag() > 5. ) { continue; }
                        if ( Z[i] < Z_min  ||  Z[i] > Z_max ) { continue; }
                        if ( p_pi[i].Theta()*rad_to_deg < 5. || p_pi[i].Theta()*rad_to_deg > 35 ){ continue; }
                	
			if( Z[i] > lead_Z_no_acc_temp ){
				lead_Z_no_acc_temp = Z[i];
				lead_idx_no_acc_temp = i;
			}	
        		

			isGoodPion_no_acc_event[i] = true;;

			double theta = p_pi[i].Theta()*rad_to_deg;

			double acc_map_pip_min = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/p_pi[i].Vect().Mag();                      
                        double acc_map_pim_min = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/p_pi[i].Vect().Mag();

			//double acc_map_pip_max = max[sector_i-1][0]->Eval(p_pi[i].Vect().Mag() );
			//double acc_map_pim_max = max[sector_i-1][1]->Eval(p_pi[i].Vect().Mag() );
			
			//double acc_map_pip_min = min[sector_i-1][0]->Eval(p_pi[i].Vect().Mag() );
			//double acc_map_pim_min = min[sector_i-1][1]->Eval(p_pi[i].Vect().Mag() );
	//	     	if ( p_pi[i].Theta()*rad_to_deg > acc_map_pip && p_pi[i].Theta()*rad_to_deg > acc_map_pim ){ 
			
			//if ( theta < acc_map_pip_max && theta > acc_map_pim_min && theta > acc_map_pip_min && theta > acc_map_pim_min ){
			if ( theta > acc_map_pip_min && theta > acc_map_pim_min ){
				if( Z[i] > lead_Z_temp ){
					lead_Z_temp = Z[i];
					lead_idx_temp = i;
				}	
		
				isGoodPion_event[i] = true;
			}

			double phi = p_pi[i].Phi()*rad_to_deg;
			double p = p_pi[i].Vect().Mag();
			//if( acceptance_match_3d( phi, theta, p, 0 ) && acceptance_match_3d( phi, theta, p, 1) ){
			if( acceptance_match_3d_cont( phi, theta, p, match3d ) ){
				if( Z[i] > lead_Z_3d_temp ){
					lead_Z_3d_temp = Z[i];
					lead_idx_3d_temp = i;
				}	
		
				isGoodPion_3d_event[i] = true;
			}
			

		}
		
		if(lead_idx_no_acc_temp > -1){
			good_events->Enter(event_count);
	
			isGoodPion_no_acc_temp.push_back(isGoodPion_no_acc_event);
			isGoodPion_temp.push_back(isGoodPion_event);
			isGoodPion_3d_temp.push_back(isGoodPion_3d_event);

			leadIdx_no_acc_vec.push_back( lead_idx_no_acc_temp );
			leadIdx_vec.push_back( lead_idx_temp );
			leadIdx_3d_vec.push_back( lead_idx_3d_temp );
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
        int nLeadIdx_3d;

 	std::vector<bool> isGoodPion;
 	std::vector<bool> isGoodPion_no_acc;
 	std::vector<bool> isGoodPion_3d;

	cout<<"Declared new variables\n";

	TTree * outTree = old_tree_upd->CloneTree(0);

	outTree->Print();        

	TBranch * branch = outTree->Branch("nLeadIdx", &nLeadIdx, "nLeadIdx/I" );
        TBranch * branch_2 = outTree->Branch("isGoodPion", &isGoodPion );
        
	TBranch * branch_3 = outTree->Branch("nLeadIdx_no_acc", &nLeadIdx_no_acc, "nLeadIdx_no_acc/I" );
        TBranch * branch_4 = outTree->Branch("isGoodPion_no_acc", &isGoodPion_no_acc);
        
	TBranch * branch_5 = outTree->Branch("nLeadIdx_3d", &nLeadIdx_3d, "nLeadIdx_3d/I" );
        TBranch * branch_6 = outTree->Branch("isGoodPion_3d", &isGoodPion_3d);

        outTree->Print();

        int nEvents = good_events->GetN();

	cout<<"Starting event loop\n";

        for(int ev = 0; ev < nEvents; ev++){

                if( ev %10000 == 0 ){cout <<ev<<" / "<<nEvents<<std::endl;}
				
                int entry  = good_events->GetEntry(ev);
                old_tree_upd->GetEntry(entry);
                
		nLeadIdx = leadIdx_vec[ev];
                nLeadIdx_no_acc = leadIdx_no_acc_vec[ev];
                nLeadIdx_3d = leadIdx_3d_vec[ev];

		isGoodPion_no_acc = isGoodPion_no_acc_temp[ev]; 
		isGoodPion_3d = isGoodPion_3d_temp[ev]; 
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

bool acceptance_match_3d( double phi_part, double theta, double p, int charge){
	
	//set momentum bins
	int this_bin_p;
	for( int i = 0; i < 4; i++ ){
		if( p > p_bin_edges[i] && p < p_bin_edges[i+1]){ this_bin_p = i; }
	}

	bool passCut = false;
	
	for( int sector = 1; sector <=6; sector++ ){
		double phi = phi_part;
		if( sector ==4 && phi < 100. ){ phi += 360; }

		//Get parameters from constants
		double theta_0 = phi_theta_bowl_theta_min[sector - 1][this_bin_p][1];
		double phi_0 = phi_theta_bowl_phi0[sector-1][this_bin_p][charge];

		//compute cut value
		double theta_min = theta_0 + pow( (phi-phi_0), 2 )/( theta_bowl_width - pow( ( phi - phi_0 ), 2 ) );

		//If phi is outside bowl, set theta_min = theta_max
		if( theta_bowl_width - pow( (phi - phi_0), 2 ) < 0 ){
			theta_min = 35.;
		}

		if( theta > theta_min ){ passCut = true; }		
	}
	return passCut;
}

bool acceptance_match_3d_cont( double phi_part, double theta, double p, TF1 * fitFuncs[6][3]){
	
	//set momentum bins
	//int this_bin_p;
	//for( int i = 0; i < 4; i++ ){
	//	if( p > p_bin_edges[i] && p < p_bin_edges[i+1]){ this_bin_p = i; }
	//}

	bool passCut[2] = {false, false};
	for( int charge = 1; charge <= 2; charge++){
		for( int sector = 0; sector < 6; sector++ ){
			double phi = phi_part;
			if( sector ==3 && phi < 100. ){ phi += 360; }

			//Get parameters from constants
			double theta_0 = fitFuncs[sector][0]->Eval(p);//phi_theta_bowl_theta_min[sector - 1][this_bin_p][1];
			double phi_0 = fitFuncs[sector][charge]->Eval(p);//phi_theta_bowl_phi0[sector-1][this_bin_p][charge];

			//compute cut value
			double theta_min = theta_0 + pow( (phi-phi_0), 2 )/( theta_bowl_width - pow( ( phi - phi_0 ), 2 ) );

			//If phi is outside bowl, set theta_min = theta_max
			if( theta_bowl_width - pow( (phi - phi_0), 2 ) < 0 ){
				theta_min = 35.;
			}
			//cout<<"theta min cont"<<theta_min<<endl;
			if( theta > theta_min ){ passCut[charge-1] = true; }		
		}
	}
	return (passCut[0] && passCut[1]);
}
