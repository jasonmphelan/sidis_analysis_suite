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

const double m_e = 0.0005;
const double m_p = 0.9383;
const double m_pi = 0.1296;

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

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<double> Ebeam_ptr(reader_rec, "Ebeam");
        TTreeReaderValue<TLorentzVector> p_e(reader_rec, "p_e");
        TTreeReaderValue<TVector3> v_e(reader_rec, "v_e");
        TTreeReaderValue<TLorentzVector> q(reader_rec, "q");

        TTreeReaderArray<TLorentzVector> p_pi(reader_rec, "p_pi");
        TTreeReaderArray<TLorentzVector> pi_q(reader_rec, "pi_q");
        TTreeReaderArray<TVector3> v_pi(reader_rec, "v_pi");
	TTreeReaderArray<int> charge(reader_rec, "charge");

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
		
		double deltaPhi_lab = 2*TMath::Pi()*(gen.Rndm());
		double deltaPhi_q = 2*TMath::Pi()*(gen.Rndm()); 
		
		TVector3 e_mom = p_e->Vect();
		TVector3 pi_mom[2];

		for( int i = 0; i < 2; i++ ){	
			//rotate pion
			TVector3 pi_q_mom = pi_q[i].Vect();
			pi_q_mom.RotateZ( deltaPhi_q );
			
			pi_mom[i] = rotate_to_beam_frame( *q, *p_e, pi_q[i] );
			pi_mom[i].RotateZ( deltaPhi_lab );

		}	

		e_mom.RotateZ(deltaPhi_lab);

	}
	
}


TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TLorentzVector pi_q ){
	TVector3 v = pi_q.Vect();
    	
	v.RotateZ( p_e.Phi() );
    	v.RotateY( q.Theta() );

    	v.RotateZ( q.Phi()  );

	return v;
}




