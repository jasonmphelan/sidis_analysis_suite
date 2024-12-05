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

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q );

int main( int argc, char** argv){

	if( argc <4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File] [Rel Uncertainty Level (%)] \n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
	double err_level = atof( argv[3] );

        TFile * file_rec = new TFile(in_name);
	TFile * outFile = new TFile(out_name, "RECREATE");
	TRandom3 gen;

	analyzer anal(0, -1);
	anal.loadAcceptanceMap( (TString)_DATA + "/acceptanceMap.root");

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
	double corr_err[2];

	outTree->Branch("beam", &beam_out);
	outTree->Branch("e", &e_out);
	outTree->Branch("pi", &pi_out);
	outTree->Branch("isGoodPion", &goodPiOut);
	outTree->Branch("Mx_2pi", &Mx_2pi_out);
	outTree->Branch("M_rho", &M_rho_out);

	outTree->Branch("rhoWeight", rhoWeight);
	outTree->Branch("rhoErr", &corr_err);

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();
		if(event_count%1000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		//initialize event variables
		bool detectedPion[2] = {true, true};
		bool goodPion[2] = {true, true};
		double onePiEvents[2] = {0};
		double twoPiEvents = 0;
		int trials = 0;

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


		//Restrict ROI
		if( !isGoodPion_acc[0] && !isGoodPion_acc[1] ){continue;} //If event wouldn't be in our final sample, continue 
		if( *Mx_2pi < 0 || *Mx_2pi > 1.7 ){continue;}
		if( *M_rho < 0 ){continue;}

		while( (corr_err[0] > err_level/100. || corr_err[1] > err_level/100.) || trials < 500){
			trials++;
			if( trials > 10000 && ( corr_err[0] > 900 || corr_err[1] > 900 ) ){ break; }
			detectedPion[0] = false;
		       	detectedPion[1] = false;
			goodPion[0] = false;
			goodPion[1] = false;

			double deltaPhi_lab = 2*TMath::Pi()*(gen.Rndm());
			double deltaPhi_q = 2*TMath::Pi()*(gen.Rndm()); 
			
			TVector3 e_mom = e->get3Momentum();
			TVector3 pi_mom[2];
		
			e_mom.RotateZ(deltaPhi_lab);
	

			bool e_acc = anal.checkAcceptance( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) ;
			if( !e_acc ){continue;}
			for( int i = 0; i < 2; i++ ){
				//rotate pion
				TVector3 pi_q_mom = pi[i].getPi_q().Vect();
				pi_q_mom.RotateZ( deltaPhi_q );
			
				pi_mom[i] = rotate_to_beam_frame( e->getQ(), e->get4Momentum(), pi[i].getPi_q() );
				pi_mom[i].RotateZ( deltaPhi_lab );
				

				//Check Pion acceptance
				int piType;
				if( pi[i].getCharge() > 0 ){ piType = 1; }
				else{ piType = 2; }
			
				int new_sec = anal.checkAcceptance( pi_mom[i].Mag(), rad_to_deg*pi_mom[i].Phi(), rad_to_deg*pi_mom[i].Theta(), piType ) ;
				if( new_sec > -1 ) detectedPion[i] = true;

				if( new_sec > -1 &&  isGoodPion[i] ){ 
					goodPion[i] =  anal.acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
				}
					
			}

			if( detectedPion[0] == true && goodPion[0] == true && detectedPion[1] == false ) {
				onePiEvents[ (int) ( pi[0].getCharge() < 0 )]++;
			}				
			if ( detectedPion[1] == true && goodPion[1] == true && detectedPion[0] == false ) {
				onePiEvents[ (int) ( pi[1].getCharge() < 0 )]++;
			}
			else if( ( detectedPion[0] == true && detectedPion[1] == true ) &&
				( goodPion[0] == true || goodPion[1] == true  ) ){
					twoPiEvents++;
			}

			//check uncertainty
			if( onePiEvents[0] != 0 && twoPiEvents != 0 ){
				corr_err[0] =  sqrt( 1./onePiEvents[0] + 1./twoPiEvents );
			}
			if( onePiEvents[1] != 0 && twoPiEvents != 0 ){
				corr_err[1] =  sqrt( 1./onePiEvents[1] + 1./twoPiEvents );
			}
			//cout<<"CURRENT UNCERTAINTY : "<<corr_err<<endl;	
	
		}	

		//Check for weird guys
		for( int i = 0; i < 2; i++ ){
			if( onePiEvents[i] == 0 || twoPiEvents == 0 ){
				rhoWeight[i] = 1.;
			}
			else{
				rhoWeight[i] =  1. + onePiEvents[i]/twoPiEvents ;
			}
		}
		outTree->Fill();

		/*
		cout<<"TRIALS FOR EVENT "<<event_count<<" : "<<trials<<endl;
 		cout<<"WEIGHT UNCERTAINTY : "<<corr_err<<endl;
		cout<<"RHO EVENT WEIGHT :"<<1+ onePiEvents/twoPiEvents<<endl;
		cout<<"\nElectron : ";
		e->get3Momentum().Print();
		cout<<"\nPi 1 : ";
		pi[0].get3Momentum().Print();
		cout<<"\nPi 2 : ";
		pi[1].get3Momentum().Print();
		cout<<"*************************************************************"<<endl;
		*/
	}
	
        cout<<"Writing to file\n";
        outFile->cd();
        outTree->Write();
        outFile->Close();
	
}


TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q ){
	TVector3 v = pi_q.Vect();
    	
	v.RotateZ( p_e.Phi() );
    	v.RotateY( q.Theta() );

    	v.RotateZ( q.Phi()  );

	return v;
}




