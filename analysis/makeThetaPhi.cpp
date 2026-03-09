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
using namespace cutVals; 

int main( int argc, char** argv){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input Path] [Output File] [# of input files] [File Type] [Beam Energy]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\nnFiles " << atoi(argv[3]) << "\n";

		TString in_name = argv[1];
       	TString out_name = argv[2];
       	int runType = atoi(argv[3]);
       	double EBeam = atof(argv[4]);

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));

	//anal.loadAcceptanceMap( (TString)_DATA + Form("/acceptance_map/acceptanceMap_%.1f.root", EBeam));
	//anal.loadMatchingFunctions();

	

	
	TChain * chain = new TChain("ePi");
	chain->Add(in_name);

	int nBins = 3;

	TH2F * hThetaPhi[2][3];

		for( int bin = 0; bin < nBins; bin++ ){
			double elMax = 35;
			double elMin = 5;
			double pimMax = 45;
			
			double p_pi = 2 + bin;
			if( p_pi>3 )pimMax = 35; // p>3
			
			hThetaPhi[0][bin] = new TH2F( Form("hThetaPhi_sec_bin_%i_pip", bin), "", 500, -250, 250, 270, 0, 35 );
			hThetaPhi[1][bin] = new TH2F( Form("hThetaPhi_sec_bin_%i_pim", bin), "", 500, -250, 250, 270, 0, pimMax );
			
		}
	

	//Load input tree
    TTreeReader reader_rec( chain );
	TTreeReaderValue<electron> e(reader_rec, "e");
    TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");
	
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
		//if( phi_e < 0 ){ phi_e += 360.; }

		int this_bin_e = (int)( ( (p_e -3.)/(8. - 3.) )*nBins);
		//if( p_e > 8 || e->getSector() == 0){continue;}

		//if( anal.checkAcceptance( e->get3Momentum().Mag(),
		//	rad_to_deg*e->get3Momentum().Phi(), 
		//	rad_to_deg*e->get3Momentum().Theta(), 0 ) < 0 ){continue;}

		//if( pass_e_fid == true ){
	

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			int piCharge = (int)( pi[i].getCharge() < 0 );
			if( !isGoodPion_no_acc[i]) continue;
			
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()*rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
			if(anal.applyAcceptanceMap( e->get3Momentum().Mag(), rad_to_deg*e->get3Momentum().Phi(), rad_to_deg*e->get3Momentum().Theta(), 0 ) <0) continue;
			if(anal.applyAcceptanceMap( p_pi, rad_to_deg*(pi)[i].get3Momentum().Phi(), rad_to_deg*(pi)[i].get3Momentum().Theta(), piCharge+1 ) < 0 ) continue;
			
			if( (pi[i].getSector() == 3 || (pi[i].getSector() == 4 && pi[i].get3Momentum().Mag() > 1) ) && phi_pi < 0. ){ phi_pi += 360; }
			
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			
			if( p_pi > 5 || pi[i].getSector() == 0 ){continue;}


			//if( anal.checkAcceptance( pi[i].get3Momentum().Mag(), 
			//	rad_to_deg*pi[i].get3Momentum().Phi(), 
			//	rad_to_deg*pi[i].get3Momentum().Theta(), 
			//	(int)( pi[i].getCharge() < 0 ) + 1 ) <0 ) { continue; }
			
			int this_bin_pi = -1;
			if(p_pi > 1.75 && p_pi < 2.25 )this_bin_pi=0;
			if(p_pi > 2.75 && p_pi < 3.25 )this_bin_pi=1;
			if(p_pi > 3.75 && p_pi < 4.25 )this_bin_pi=2;
			if( this_bin_pi < 0 ) continue;

			if (this_bin_pi < 0 || this_bin_pi > 2) {continue;}
			hThetaPhi[chargeIdx][this_bin_pi]->Fill(phi_pi, theta_pi);
			
		}
	}

	cout<<"FINISHED LOOP\n";
        TFile * outFile = new TFile(out_name, "RECREATE");
	outFile->cd();	
	
		for( int bin = 0; bin < nBins; bin++ ){
			hThetaPhi[0][bin]->Write();
			hThetaPhi[1][bin]->Write(); 
		}
	
	
	outFile->Close();
}

