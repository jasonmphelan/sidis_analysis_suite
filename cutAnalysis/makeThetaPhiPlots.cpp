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

	analyzer anal(0, -1);
	//anal.loadAcceptanceMap( (TString)_DATA + Form("/acceptance_map/acceptanceMap_%.1f.root", EBeam));
	//anal.loadMatchingFunctions();

	cout<<"FILE MADE\n";
	reader skimReader;
	skimReader.setNumFiles( nFiles);
	skimReader.setRunType( runType );
	skimReader.setEnergy( EBeam );

	cout<<"GETTING FILES \n";
	TChain * chain = new TChain("ePi");
	skimReader.getRunSkimsByName(chain, in_name);

	int nBins = 40;

	TH2F * hThetaPhi[3][6][nBins];


	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < nBins; bin++ ){
			double elMax = 35;
			double elMin = 5;
			double pimMax = 45;
			double p_e = 3 + 5.*bin/(double)nBins;
			double p_pi = 0 + 5.*bin/(double)nBins;
			if( p_pi>3 )pimMax = 35; // p>3
			if( p_e >= 3.75 && p_e < 4.5 ) elMax = 30;
			if( p_e >= 5.5 && p_e < 6.75) elMax = 25;
			if( p_e >= 6.75 ) elMax = 15;
			hThetaPhi[1][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_pip", sec, bin), "", 500, -250, 250, 270, 0, 35 );
			hThetaPhi[2][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_pim", sec, bin), "", 500, -250, 250, 270, 0, pimMax );
			hThetaPhi[0][sec][bin] = new TH2F( Form("hThetaPhi_sec_%i_bin_%i_e", sec, bin), "", 500, -250, 250, 270, 0, elMax );
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

		bool pass_e_fid = true;
		//for (int regionIdx=0; regionIdx<3; regionIdx++) {
		//	if( e->getEdge(regionIdx) < 1.5*e_fid_cuts[regionIdx] ){ pass_e_fid = false; }
		//}	
		


		double p_e = e->get3Momentum().Mag();
		double theta_e = e->get3Momentum().Theta()*rad_to_deg;
		double phi_e = e->get3Momentum().Phi()*rad_to_deg;
		if( e->getSector() == 4 && phi_e < 100. ){ phi_e += 360; }
		//if( phi_e < 0 ){ phi_e += 360.; }

		int this_bin_e = (int)( ( (p_e -3.)/(8. - 3.) )*nBins);
		//if( p_e > 8 || e->getSector() == 0){continue;}

		//if( anal.checkAcceptance( e->get3Momentum().Mag(),
		//	rad_to_deg*e->get3Momentum().Phi(), 
		//	rad_to_deg*e->get3Momentum().Theta(), 0 ) < 0 ){continue;}

		//if( pass_e_fid == true ){
		if( p_e > 3 && p_e < 8 && e->getSector() != 0){
			hThetaPhi[0][e->getSector()-1][this_bin_e]->Fill(phi_e, theta_e);
		}
		

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			bool pass_fiducials = true;
			int piCharge = (int)( pi[i].getCharge() < 0 );

			//for (int regionIdx=0; regionIdx<3; regionIdx++) {
			//	if( pi[i].getEdge(regionIdx) < 1.5*pi_fid_cuts[piCharge][regionIdx] ) pass_fiducials =false;
			//}
		
			//if( pass_fiducials == false ) continue;

			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()*rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
			if( (pi[i].getSector() == 3 || (pi[i].getSector() == 4 && pi[i].get3Momentum().Mag() > 1) ) && phi_pi < 0. ){ phi_pi += 360; }
			
			int chargeIdx = (int)( pi[i].getCharge() < 0 ) + 1;
			
			if( p_pi > 5 || pi[i].getSector() == 0 ){continue;}

			//if( anal.checkAcceptance( pi[i].get3Momentum().Mag(), 
			//	rad_to_deg*pi[i].get3Momentum().Phi(), 
			//	rad_to_deg*pi[i].get3Momentum().Theta(), 
			//	(int)( pi[i].getCharge() < 0 ) + 1 ) <0 ) { continue; }
			
			int this_bin_pi = (int)( ( (p_pi - 0.)/(5.- 0.) )*nBins);
			
			hThetaPhi[chargeIdx][pi[i].getSector()-1][this_bin_pi]->Fill(phi_pi, theta_pi);
			
		}
	}

	cout<<"FINISHED LOOP\n";
        TFile * outFile = new TFile(out_name, "RECREATE");
	outFile->cd();	
	
	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < nBins; bin++ ){
			hThetaPhi[0][sec][bin]->Write();
		}
		for( int bin = 0; bin < nBins; bin++ ){
			hThetaPhi[1][sec][bin]->Write();
			hThetaPhi[2][sec][bin]->Write(); 
		}
	}
	
	outFile->Close();
}

