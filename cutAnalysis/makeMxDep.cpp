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

	if( argc < 5 ){
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

	cout<<"Read inputs\n";

	reader skimReader;
	skimReader.setNumFiles( nFiles);
	skimReader.setRunType( runType );
	skimReader.setEnergy( EBeam );

	cout<<"Set Reader\n";

	TChain * chain = new TChain("ePi");
	skimReader.getRunSkimsByName(chain, in_name);

	cout<<"Loaded skims\n";

	TString data_type[2] = {"pip", "pim"};
	TH1F* h_W[2];    
	TH1F* h_Xb[2];    
	TH1F* h_Q2[2];     
	TH1F* h_eta[2];     
	TH1F* h_y[2];        
	TH1F* h_omega[2];     
	TH1F* h_Pt_e[2];       
	TH1F* h_theta_e[2];     
	TH1F* h_Vz_e[2];         
	TH1F* h_p_e[2];           
	TH1F* h_phi_e[2];        

	TH1F* h_Z[2];             
	TH1F* h_p_pi[2];          
	TH1F* h_Vz_pi[2];         
	TH1F* h_delta_phi_pi[2];  
	TH1F* h_theta_pi[2];     
	TH1F* h_phi_pi[2];        
	TH1F* h_Pt_pi[2];
	TH1F* h_Mx_Q2[2][4]; 
	TH1F* h_Mx_Z[2][4]; 
	

	double q2_bin_edge[5] = {2, 3, 4, 6, 8};
	double z_bin_edge[5] = {.3, .4, .5, .6, .8};

	for( int i = 0; i < 2; i++ ){
		int j = 0;
		int k = 0;
		h_Q2[i]            = new TH1F("hQ2_"+data_type[i], Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 200, 2, 8);
	
		h_Z[i]             = new TH1F("hZ_"+data_type[i], Form("Z_%i_%i;Z;Counts [a.u.]", j, k), 100, 0.3, .8);
		
		for( int j = 0; j < 4; j++ ){
			h_Mx_Q2[i][j]            = new TH1F("hMx_Q2"+data_type[i] + Form("_%i", j), Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 300, 0, 3);
			h_Mx_Z[i][j]            = new TH1F("hMx_Z"+data_type[i] + Form("_%i", j), Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 300, 0, 3);
		}
	}

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);

	//Load input tree
        TTreeReader reader_rec( chain );
	TTreeReaderValue<electron> e(reader_rec, "e");
        TTreeReaderArray<pion> pi(reader_rec, "pi");

	cout<<"BEGIN EVENT LOOP\n";
	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}


		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			int idx = (int)( pi[i].getCharge() < 0 );
	


			if( sqrt(e->getW2()) < 2.5  ){continue;}
			
			if( e->getQ2() < 2 ){ continue; }

			if( e->getY() > 0.75 ){continue;}

			if( pi[i].get3Momentum().Mag() < 1.25 || pi[i].get3Momentum().Mag() > 5 ){continue;}
			if( e->get3Momentum().Theta()*rad_to_deg < 5 || e->get3Momentum().Theta()*rad_to_deg > 35 ){continue;}
		//	if( pi[i].get3Momentum().Theta()*rad_to_deg < 5 || pi[i].get3Momentum().Theta()*rad_to_deg > 35 ){continue;}
			if( pi[i].getZ() < 0.3 ){ continue; }	
			
			int q2 = -1;
			int z = -1;

			for( int bin = 0; bin < 4; bin++ ){
				if( e->getQ2() > q2_bin_edge[bin] && e->getQ2() < q2_bin_edge[bin+1] ) q2 = bin;
				if( pi[i].getZ() > z_bin_edge[bin] && pi[i].getZ() < z_bin_edge[bin+1] ) z = bin;
			}
			if( q2 < 0 || z < 0 ){ continue; }

			h_Mx_Q2[idx][q2]->Fill( pi[i].getMx() );
			h_Mx_Z[idx][z]->Fill( pi[i].getMx() );
			h_Q2[idx]->Fill( e->getQ2() );
			h_Z[idx]->Fill( pi[i].getZ() );
		}
	}
	
	cout<<"Completed good event list... \n";

        cout<<"Writing to file\n";
        
	TFile * outFile = new TFile( out_name, "RECREATE");
	outFile->cd();
	for( int i = 0; i < 2; i++ ){
		h_Q2[i]->Write();            
		h_Z[i]->Write();             
	
		for( int j = 0; j < 4; j++ ){
			h_Mx_Q2[i][j]->Write();            
			h_Mx_Z[i][j]->Write();            
		}
	}

        outFile->Close();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout<<"Done. Elapsed Time : "<<elapsed.count()<<std::endl;

	return 0;
}
