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

	cout<<"FILE MADE\n";
	reader skimReader;
	skimReader.setNumFiles( nFiles);
	skimReader.setRunType( runType );
	skimReader.setEnergy( EBeam );

	cout<<"GETTING FILES \n";
	TChain * chain = new TChain("ePi");
	skimReader.getRunSkimsByName(chain, in_name);


	TH2F * hThetaPhi[3][6][10];

	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < 5; bin++ ){
			hThetaPhi[1][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_pip", sec, bin), "", 880, -220, 220, 80, 0, 40 );
			hThetaPhi[2][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_pim", sec, bin), "", 880, -220, 220, 80, 0, 40 );
		}
		for( int bin = 0; bin < 10; bin++ ){
			hThetaPhi[0][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_e", sec, bin), "", 880, -220, 220, 80, 0, 40 );
		}
	}


	//Load input tree
        TTreeReader reader_rec( chain );
	TTreeReaderValue<electron> e(reader_rec, "e");
        TTreeReaderArray<pion> pi(reader_rec, "pi");
	
	int event_total = reader_rec.GetEntries();
	cout<<"BEGIN EVENT LOOP\n";	
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		double p_e = e->get3Momentum().Mag();
		double theta_e = e->get3Momentum().Theta()*rad_to_deg;
		double phi_e = e->get3Momentum().Phi()*rad_to_deg;
		//if( phi_e < 0 ){ phi_e += 360.; }

		int this_bin_e = (int)( ( (p_e -0.)/(10. - 0.) )*10.);
		if( p_e > 10 ){continue;}


		hThetaPhi[0][e->getSector()-1][this_bin_e]->Fill(phi_e, theta_e);

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()*rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
			//if( phi_pi < 0 ){ phi_pi += 360.; }
			
			int chargeIdx = (int)( pi[i].getCharge() < 0 ) + 1;

			
			int this_bin_pi = (int)( ( (p_pi - 0.)/(5.- 0.) )*5.);
			if( p_pi > 5 ){continue;}
			hThetaPhi[chargeIdx][pi[i].getSector()-1][this_bin_pi]->Fill(phi_pi, theta_pi);
			
		}
	}

	cout<<"FINISHED LOOP\n";
        TFile * outFile = new TFile(out_name, "RECREATE");
	outFile->cd();	
	
	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < 10; bin++ ){
			hThetaPhi[0][sec][bin]->Write();
		}
		for( int bin = 0; bin < 5; bin++ ){
			hThetaPhi[1][sec][bin]->Write();
			hThetaPhi[2][sec][bin]->Write(); 
		}
	}
	
	outFile->Close();
}

