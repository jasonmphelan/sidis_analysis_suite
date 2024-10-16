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


int main( int argc, char** argv){

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Pi File] [K File]  [Output File] [Matching]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2]  << "\n";

	TString inName_rec = argv[1];
	TString inName_gen = argv[2];
	TString outName = argv[3];

	TFile * outFile = new TFile( outName, "RECREATE");

	TFile * inFile_rec = new TFile( inName_rec );
	TTree * piChain = (TTree *)inFile_rec->Get("ePi");
	
	TFile * inFile_gen = new TFile( inName_gen );
	TTree * kChain = (TTree *)inFile_gen->Get("ePi");

	TH1F * kHists[bins_xB][bins_Q2][2];
	TH1F * piHists[bins_xB][bins_Q2][2];
	

	for( int i = 0; i < bins_xB; i++ ){
		for( int j = 0; j < bins_Q2; j++ ){
		
			kHists[i][j][0] = new TH1F( (TString)"kHist_P_"+Form("_%i_%i", i, j), Form("kHist_P_%i_%i", i, j), 10, .3, .8); 
			piHists[i][j][0] = new TH1F((TString)"piHist_P_"+Form("_%i_%i", i, j), Form("piHist_P_%i_%i", i, j), 10, .3, .8); 
			
			kHists[i][j][1] = new TH1F( (TString)"kHist_M_"+Form("_%i_%i", i, j), Form("kHist_M_%i_%i", i, j), 10, .3, .8); 
			piHists[i][j][1] = new TH1F( (TString)"piHist_M_"+Form("_%i_%i", i, j), Form("piHist_M_%i_%i", i, j), 10, .3, .8); 
			
		}
	}

        TTreeReader reader_pi(piChain);

        TTreeReaderValue<electron> e_pi_ptr(reader_pi, "e");
        TTreeReaderArray<pion> pi(reader_pi, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_pi, "isGoodPion");

	int event_total = piChain->GetEntries();

	while (reader_pi.Next()) {
                int event_count = reader_pi.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

                double Q2 = e_pi_ptr->getQ2();
                double xB = e_pi_ptr->getXb();
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
                int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB);

			
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			if( !isGoodPion[i] ){continue;}

			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			piHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi[i].getZ() );
			
		}
	}
				
        TTreeReader reader_k(kChain);

        TTreeReaderValue<electron> e_k_ptr(reader_k, "e");
        TTreeReaderArray<pion> kaon(reader_k, "pi");
	TTreeReaderArray<bool> isGoodKaon(reader_k, "isGoodPion");

	event_total = kChain->GetEntries();

	while (reader_k.Next()) {
                int event_count = reader_k.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

                double Q2 = e_k_ptr->getQ2();
                double xB = e_k_ptr->getXb();
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
                int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB);

			
		for( int i = 0; i < (int) ( kaon.end() - kaon.begin() ); i++ ){
			if( !isGoodKaon[i] ){continue;}

			int chargeIdx = (int)( kaon[i].getCharge() < 0 );
			kHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( kaon[i].getZ() );
			
		}
	}


	
	

	outFile->cd();


	for( int i = 0; i < bins_xB; i++ ){
		for( int j = 0; j < bins_Q2; j++ ){
			
			kHists[i][j][0]->Divide(piHists[i][j][0]);	
			kHists[i][j][1]->Divide(piHists[i][j][1]);	
					
			kHists[i][j][0]->Write();
			kHists[i][j][1]->Write();

			
		}
	}
	outFile->Close();
	//file_1->Close();
	//file_2->Close();
	
	return 1;
}
