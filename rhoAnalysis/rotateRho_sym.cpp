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
#include "TRandom3.h"
#include "CLAS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "constants.h"

using namespace constants;

#define CORR_PATH _DATA
const int thread_MAX = 50;

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TVector3 pi_q );

int main( int argc, char** argv){

	if( argc < 6 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File] [beam energy] [Rel Uncertainty Level (%)] \n";
		return -1;
	}

	TString in_name = argv[1];
    TString out_name = argv[2];
	double energy = atof( argv[3] );
	double err_level = atof( argv[4] );
	int acc_match = atoi( argv[5] );
	TString acc_name = argv[6];
	int i_thread = atoi( argv[7] );

    TFile * file_rec = new TFile(in_name);
	TFile * outFile = new TFile(out_name  + Form("_chunk_%i.root", i_thread), "RECREATE");
	TRandom3 gen(0);

	analyzer anal(0, -1);
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/"+acc_name);//%.1f.root", energy));
	anal.loadMatchingFunctions("matchCut2D_map.root");//(TString)_DATA + (TString)"/acceptance_matching/matchCut2D_map.root");
	
	
	//Load input tree
    TTreeReader reader_rec("ePi", file_rec);
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion_no_acc");
	TTreeReaderArray<bool> isGoodPion_acc(reader_rec, "isGoodPion");
	TTreeReaderValue<TLorentzVector> beam(reader_rec, "beam");

	//Set output tree
	TTree * outTree = new TTree("ePi", "(e,e'pi) event  information");
	TLorentzVector beam_out;	
	electron e_out;
	std::vector<pion> pi_out;
	std::vector<bool> goodPiOut;
	
	double Mx_2pi_out;
	double M_rho_out;
	double rhoWeight[2];
	double rhoWeight_sym[2];
	double corr_err[2];
	double corr_err_sym[2];
	int trials;
	int event_num;

	outTree->Branch("beam", &beam_out);
	outTree->Branch("e", &e_out);
	outTree->Branch("pi", &pi_out);
	outTree->Branch("isGoodPion", &goodPiOut);
	outTree->Branch("Mx_2pi", &Mx_2pi_out);
	outTree->Branch("M_rho", &M_rho_out);

	outTree->Branch("rhoWeight", rhoWeight, "rhoWeight[2]/D");
	outTree->Branch("rhoWeight_sym", rhoWeight_sym, "rhoWeight_sym[2]/D");

	outTree->Branch("rhoErr", corr_err, "rhoErr[2]/D");
	outTree->Branch("rhoErr_sym", corr_err_sym, "rhoErr[2]/D");

	outTree->Branch("trials", &trials);
	outTree->Branch("event_num", &event_num);
	std::ofstream txtFile;
	txtFile.open( "../plotting/rho_data_points_one_pi.txt" );
	
	
	int event_total = reader_rec.GetEntries();
	int chunkSize  = event_total / (thread_MAX-1);
	auto start = std::chrono::high_resolution_clock::now();

	while (reader_rec.Next()) {
        int event_count = reader_rec.GetCurrentEntry();
		if( event_count < i_thread*chunkSize || event_count >= (i_thread+1)*chunkSize )continue;

		event_num = event_count;
		
		if((event_count%500) == (0)){
			cout<<"Thread "<<i_thread<<" : "<<100.*(double)(event_count - i_thread*chunkSize)/(double)chunkSize<<"%";
			auto finish = std::chrono::high_resolution_clock::now();
    		std::chrono::duration<double> elapsed = finish - start;
			cout<<" in "<<elapsed.count()/60.<<" minutes\n";
		}		
	
		//initialize event variables
		
		bool detectedPion[2] = {false, false};
		bool goodPion[2] = {false, false};
		bool accPion[2] = {false, false};

		bool detectedPion_sym[2] = {false, false};
		bool goodPion_sym[2] = {false, false};
		bool accPion_sym[2] = {false, false};

		double onePiEvents[2] = {0};
		double onePiEvents_sym[2] = {0};

		double twoPiEvents = 0;
		double twoPiEvents_sym = 0;

		double zeroPiEvents = 0;
		trials = 0;

		pi_out.clear();
		goodPiOut.clear();
		e_out.Clear();

		beam_out = (TLorentzVector) (*beam);
		e_out = (electron)(*e);
		pi_out.push_back((pion)pi[0]);
		pi_out.push_back((pion)pi[1]);
		
		goodPiOut.push_back((bool)isGoodPion[0]);
		goodPiOut.push_back((bool)isGoodPion[1]);
		
		Mx_2pi_out = (double)(*Mx_2pi);
		M_rho_out = (double)(*M_rho);
		rhoWeight[0] = 0; // (pi+, pi-)
		rhoWeight[1] = 0; // (pi+, pi-)
		corr_err[0] = 999;
		corr_err[1] = 999;

		rhoWeight_sym[0] = 0; // (pi+, pi-)
		rhoWeight_sym[1] = 0; // (pi+, pi-)
		corr_err_sym[0] = 999;
		corr_err_sym[1] = 999;

		double p_e = e->get3Momentum().Mag();
		double theta_e = e->get3Momentum().Theta()*rad_to_deg;
		double phi_e = e->get3Momentum().Phi()*rad_to_deg;

		//Look at particles in acceptance map
		if( anal.applyAcceptanceMap( p_e,
					phi_e, 
					theta_e, 0 )  < 0){continue;}
		if( anal.applyAcceptanceMap( pi[0].get3Momentum().Mag(), 
					rad_to_deg*pi[0].get3Momentum().Phi(), 
					rad_to_deg*pi[0].get3Momentum().Theta(), 
					(int)( pi[0].getCharge() < 0 ) + 1 ) < 0 ) { continue; }
		if( anal.applyAcceptanceMap( pi[1].get3Momentum().Mag(), 
					rad_to_deg*pi[1].get3Momentum().Phi(), 
					rad_to_deg*pi[1].get3Momentum().Theta(), 
					(int)( pi[1].getCharge() < 0 ) + 1 ) < 0  ) { continue; }
		//Restrict ROI			
		if( acc_match && !anal.applyAcceptanceMatching(pi[0], 2) ){continue;}
		if( acc_match && !anal.applyAcceptanceMatching(pi[1], 2) ){continue;}
		//if( !isGoodPion_acc[0] && !isGoodPion_acc[1] ){continue;} //If event wouldn't be in our final sample, continue 
		if( *Mx_2pi < 0 || *Mx_2pi > 2.5 ){continue;}
		if( *M_rho < 0 ){continue;}

		int symType = 0;

		int symSec_0 = anal.applyAcceptanceMap( pi[0].get3Momentum().Mag(), 
						rad_to_deg*pi[0].get3Momentum().Phi(), 
						rad_to_deg*pi[0].get3Momentum().Theta(), 
						(int)( pi[1].getCharge() < 0 ) + 1 ) ;
		int symSec_1 = anal.applyAcceptanceMap( pi[1].get3Momentum().Mag(), 
						rad_to_deg*pi[1].get3Momentum().Phi(), 
						rad_to_deg*pi[1].get3Momentum().Theta(), 
						(int)( pi[0].getCharge() < 0 ) + 1 );
		

		if ( symSec_0 >= 0 && symSec_1 >= 0)symType = 2;
		else if (symSec_0 >= 0 && symSec_1 < 0) symType = 0;
		else if (symSec_1 >= 0 && symSec_0 < 0) symType = 1;
		else{ symType = -1; }
		//Find electron rotation angles that satisfy acceptance map for fixed p, theta

		std::vector<std::pair<double, double>> phi_ranges;
		std::vector<double> cumulative_phi;
		double total_phi = 0;

		for( int sec=0 ; sec < 6; sec++ ){
			double theta_min = anal.mapFunc(0, sec, 4, p_e);
			double theta_max = anal.mapFunc(0, sec, 5, p_e);

			if(theta_e > theta_min && theta_e < theta_max){
				phi_ranges.push_back({anal.returnElPhiMin(sec, p_e, theta_e), anal.returnElPhiMax(sec, p_e, theta_e)});
				///cout<<anal.returnElPhiMin(sec, p_e, theta_e)<< " and "<<anal.returnElPhiMax(sec, p_e, theta_e)<<std::endl;
				total_phi += (anal.returnElPhiMax(sec, p_e, theta_e) - anal.returnElPhiMin(sec, p_e, theta_e));
				cumulative_phi.push_back(total_phi);
				
			}
		}


		//Keep rotating while (err > input) for each pion.  Also set a mininum number of trials
		//f( event_count >500 )break;
		while( (isGoodPion[0] && (corr_err[0] > err_level/100. || corr_err_sym[0] > err_level/100.)) || (isGoodPion[1] && ( corr_err[1] > err_level/100. || corr_err_sym[1] > err_level/100. )) || trials < 1000){
			trials++;

			if ( trials >500000  )break;

			//if( trials > 50000 && 
			//	( ( isGoodPion[0] && corr_err[0] > 900) || (isGoodPion[1] && corr_err[1] > 900 )
			 //) ){ break; }
			
			double symType = 0;

			detectedPion[0] = false;
		    detectedPion[1] = false;
			goodPion[0] = false;
			goodPion[1] = false;
			accPion[0] = false;
			accPion[1] = false;

			detectedPion_sym[0] = false;
		    detectedPion_sym[1] = false;
			goodPion_sym[0] = false;
			goodPion_sym[1] = false;
			accPion_sym[0] = false;
			accPion_sym[1] = false;
		
		

			//initialize rotation angles
			double newPhi_lab = 0;
			double r = gen.Uniform(0, total_phi);
			for ( size_t i = 0; i < phi_ranges.size(); i++ ){
				//cout<<cumulative_phi[i]<<std::endl;
				if( r < cumulative_phi[i] ){
					double base = (i==0)?0.0 : cumulative_phi[i-1];

					
					double offset = r-base;
					newPhi_lab = phi_ranges[i].first + offset;
					break;

				}
			}

			double deltaPhi_lab =  (newPhi_lab - phi_e)/rad_to_deg;//2*TMath::Pi()*(gen.Rndm());//(newPhi_lab - phi_e)/rad_to_deg;//
			double deltaPhi_q = 2*TMath::Pi()*(gen.Rndm()); 
			
			TVector3 e_mom = e->get3Momentum();
			TVector3 pi_mom[2];

		
		
			//rotate electron about beam axis and check acceptance
			//e_mom.RotateZ(deltaPhi_lab);
			//We can choose only angles that would already satisfy acceptance map
			//int e_acc = anal.applyAcceptanceMap( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) ;
			//cout<<e_acc<<std::endl;
			//if( e_acc < 0 ){continue;}

			for( int i = 0; i < 2; i++ ){
				//rotate pion about q and z
				TVector3 pi_q_mom = pi[i].getPi_q().Vect();
				pi_q_mom.RotateZ( deltaPhi_q );
			
				pi_mom[i] = rotate_to_beam_frame( e->getQ(), e->get4Momentum(), pi_q_mom );
				pi_mom[i].RotateZ( deltaPhi_lab );
				

				//Check Pion acceptance
				int piType = (int)(pi[i].getCharge() < 0) + 1;
			
				
				int new_sec = anal.applyAcceptanceMap( pi_mom[i].Mag(), rad_to_deg*pi_mom[i].Phi(), rad_to_deg*pi_mom[i].Theta(), piType ) ;
				
			
		

				if( new_sec > -1 ) {
					if( acc_match ) detectedPion[i] = anal.acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
					else detectedPion[i] = true;
				}
				
				if( new_sec > -1 &&  isGoodPion[i] ){ 					
					accPion[i] = true;
					goodPion[i] =  anal.acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
				}


				if( symType == 0 || symType == 1 )continue;
				int new_sec_sym = anal.applyAcceptanceMap( pi_mom[i].Mag(), rad_to_deg*pi_mom[i].Phi(), rad_to_deg*pi_mom[i].Theta(), !((bool)(piType - 1) ) + 1) ;			
				//cout<<"Not one pi";
				if( new_sec_sym > -1 ) {
					if( acc_match ) detectedPion_sym[i] = anal.acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
					else detectedPion_sym[i] = true;
				}
				if( new_sec_sym > -1 &&  isGoodPion[i] ){ 					
					accPion_sym[i] = true;
					goodPion_sym[i] =  anal.acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
				}
			}

			if( detectedPion[0] == true && goodPion[0] == true && detectedPion[1] == false ) {
				onePiEvents[ (int) ( pi[0].getCharge() < 0 )]++;
			}				
			if ( detectedPion[1] == true && goodPion[1] == true && detectedPion[0] == false ) {
				onePiEvents[ (int) ( pi[1].getCharge() < 0 )]++;
			}
			if( ( detectedPion[0] == true && detectedPion[1] == true ) &&
				( goodPion[0] == true || goodPion[1] == true  ) ){			
					twoPiEvents++;
			}

			if( detectedPion_sym[0] == true && goodPion_sym[0] == true && detectedPion_sym[1] == false ) {
				onePiEvents_sym[ (int) ( pi[0].getCharge() < 0 )]++;
			}				
			if ( detectedPion_sym[1] == true && goodPion_sym[1] == true && detectedPion_sym[0] == false ) {
				onePiEvents_sym[ (int) ( pi[1].getCharge() < 0 )]++;
			}
			if( ( detectedPion_sym[0] == true && detectedPion_sym[1] == true ) &&
				( goodPion_sym[0] == true || goodPion_sym[1] == true  ) ){			
					twoPiEvents_sym++;
			}
			if( ( detectedPion_sym[0] == false && detectedPion_sym[1] == false ) ){			
					zeroPiEvents++;
			}



			//check uncertainty
			if( onePiEvents[0] != 0 && twoPiEvents != 0 ){
				corr_err[0] =  sqrt( 1./onePiEvents[0] + 1./twoPiEvents );
			}
			if( onePiEvents[1] != 0 && twoPiEvents != 0 ){
				corr_err[1] =  sqrt( 1./onePiEvents[1] + 1./twoPiEvents );
			}

			if( symType == 0 || symType == 1 ) {
				corr_err_sym[0] = 0;
				corr_err_sym[1] = 0;
			}
			double denom = 0;
			if( symType == 2 ){ denom = twoPiEvents_sym; }
			if( symType == -1 ){ denom = zeroPiEvents; }
			if( denom > 0 ){
				if( onePiEvents_sym[0] != 0 && twoPiEvents_sym != 0 ){
					corr_err_sym[0] =  sqrt( 1./onePiEvents_sym[0] + 1./denom );
				}
				if( onePiEvents_sym[1] != 0 && twoPiEvents_sym != 0 ){
					corr_err_sym[1] =  sqrt( 1./onePiEvents_sym[1] + 1./denom );
				}
			}
			//cout<<"CURRENT UNCERTAINTY : "<<corr_err<<endl;	
	
		}	
		if( event_count == 193)txtFile.close();


		//Check for weird guys
		for( int i = 0; i < 2; i++ ){
			if( onePiEvents[i] == 0 || twoPiEvents == 0 ){
				if( isGoodPion_acc[i] ){
					rhoWeight[i] = 1.;
				}
				else{ rhoWeight[i] = 0.;}
			}
			else if( isGoodPion_acc[i] ){
				rhoWeight[i] =  1. + onePiEvents[i]/twoPiEvents ;
			}
			else{
				rhoWeight[i] = onePiEvents[i]/twoPiEvents ;
			}
		}
		if( symType == 0 ){
			rhoWeight_sym[0] = 1.;
			rhoWeight_sym[1] = 0.;
		}
		else if( symType == 1 ){
			rhoWeight_sym[1] = 1.;
			rhoWeight_sym[0] = 0.;
		}
		else if( symType == 2 ){
			rhoWeight_sym[0] =  1. + onePiEvents_sym[0]/twoPiEvents_sym ;
			rhoWeight_sym[1] =  1. + onePiEvents_sym[1]/twoPiEvents_sym ;
		}
		else{
			rhoWeight_sym[0] =  1. + onePiEvents_sym[0]/zeroPiEvents ;
			rhoWeight_sym[1] =  1. + onePiEvents_sym[1]/zeroPiEvents ;
		}



		
		
		outTree->Fill();	
	}    
    outFile->cd();
    outTree->Write();
    outFile->Close();
		
	
}



TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TVector3 pi_q ){
	TVector3 v = pi_q;
    
	p_e.RotateZ(-q.Phi());
	p_e.RotateY(-q.Theta());

	v.RotateZ( p_e.Phi() );
    v.RotateY( q.Theta() );

    v.RotateZ( q.Phi()  );

	return v;
}




