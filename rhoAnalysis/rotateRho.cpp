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

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q );

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

        TFile * file_rec = new TFile(in_name);
	TRandom3 gen;

	CLAS acceptanceMap;
	analyzer anal(0, -1);
	anal.loadMatchingFunctions();
	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
		
		double deltaPhi_lab = 2*TMath::Pi()*(gen.Rndm());
		double deltaPhi_q = 2*TMath::Pi()*(gen.Rndm()); 
		
		TVector3 e_mom = e->get3Momentum();
		TVector3 pi_mom[2];
		
		e_mom.RotateZ(deltaPhi_lab);

		int e_acc = acceptanceMap.GetElectronAcceptance( e_mom.Theta(), e_mom.Phi(), e_mom.Mag() ) ;
		cout<<"Electron : "<<e_acc<<endl;
		

		for( int i = 0; i < 2; i++ ){	
			//rotate pion
			TVector3 pi_q_mom = pi[i].getPi_q().Vect();
			pi_q_mom.RotateZ( deltaPhi_q );
			
			pi_mom[i] = rotate_to_beam_frame( e->getQ(), e->get4Momentum(), pi[i].getPi_q() );
			pi_mom[i].RotateZ( deltaPhi_lab );
			int pi_acc = anal.acceptance_match_3d_cont( pi_mom[i].Phi()* rad_to_deg, pi_mom[i].Theta()*rad_to_deg , pi_mom[i].Mag(), (int) ( pi[i].getCharge() < 0 ));
			cout<<"Pion #"<<i+1<<" : "<<pi_acc<<endl;

		}	

		cout<<"-------------------------------------------------------------------------\n";
	}
	
}


TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q ){
	TVector3 v = pi_q.Vect();
    	
	v.RotateZ( p_e.Phi() );
    	v.RotateY( q.Theta() );

    	v.RotateZ( q.Phi()  );

	return v;
}




