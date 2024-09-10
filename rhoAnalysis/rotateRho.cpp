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
		cerr << "./code [Input File] [Output File] [Rel Uncertainty Level] \n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
	double err_level = atof( argv[3] );

        TFile * file_rec = new TFile(in_name);
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

	double onePiEvents = 0;
	double twoPiEvents = 0;
	double corr_err = 999;
	bool detectedPion[2] = {true, true};

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();
		int trials = 0;

		if(event_count%1000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
		if( *Mx_2pi > 1.2 || *Mx_2pi < .75 ){continue;}
		if( *M_rho > 1 || *M_rho < .45 ){continue;}
		while( corr_err > err_level ){
			trials++;
			detectedPion[0] = false;
		       	detectedPion[1] = false;


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
				//detectedPion[i] = anal.checkPiAcceptance( 	
				int piType;
				if( pi[i].getCharge() > 0 ){ piType = 1; }
				else{ piType = 2; }
				detectedPion[i] = anal.checkAcceptance( pi_mom[i].Mag(), rad_to_deg*pi_mom[i].Phi(), rad_to_deg*pi_mom[i].Theta(), piType ) ;
			}

			if( ( detectedPion[0] == true && isGoodPion[0] == true && detectedPion[1] == false ) ||
				( detectedPion[1] == true && isGoodPion[1] == true && detectedPion[0] == false ) ){
					onePiEvents++;

					//cout<<"FOUND 1 PI EVENT\n";
			}
			else if( ( detectedPion[0] == true && detectedPion[1] == true ) &&
				( isGoodPion[0] == true || isGoodPion[1] == true  ) ){
					twoPiEvents++;
					//cout<<"FOUND 2 PI EVENT! \n";
			}

			//check uncertainty
			if( onePiEvents != 0 && twoPiEvents != 0 ){
				corr_err =  sqrt( 1/onePiEvents + 1/twoPiEvents );
			}
			cout<<"CURRENT UNCERTAINTY : "<<corr_err<<endl;	
		}	
		cout<<"TRIALS FOR EVENT "<<event_count<<" : "<<trials<<endl;
 		cout<<"RHO EVENT WEIGHT :"<<1+ onePiEvents/twoPiEvents<<endl;

	}
	
}


TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q ){
	TVector3 v = pi_q.Vect();
    	
	v.RotateZ( p_e.Phi() );
    	v.RotateY( q.Theta() );

    	v.RotateZ( q.Phi()  );

	return v;
}




