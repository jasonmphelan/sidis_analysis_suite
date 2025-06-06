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
	TH1F* h_Mx[2]; 
	
	TH2F* hQ2_W[2];
	TH2F* hQ2_omega[2];

	TH2F* hZ_Mx[2];
	TH2F* hQ2_Mx[2];
	TH2F* hXb_Mx[2];
	TH2F* hPt_Mx[2];


	for( int i = 0; i < 2; i++ ){
		int j = 0;
		int k = 0;
		h_W[i]             = new TH1F("hW_"+data_type[i], Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k), 250, 0, 4.75);
		h_Xb[i]            = new TH1F("hXb_"+data_type[i], Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k), 50, 0, 1);
		h_Q2[i]            = new TH1F("hQ2_"+data_type[i], Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 250, 0, 10);
		h_eta[i]           = new TH1F("hEta_"+data_type[i], Form("eta_%i_%i;#eta;Counts [a.u.]", j, k), 50, 0, 4.6);
		h_y[i]             = new TH1F("hY_"+data_type[i], Form("y_%i_%i;y;Counts [a.u.]", j, k), 250, 0, 1);
		h_omega[i]         = new TH1F("hOmega_"+data_type[i], Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k), 50, 0, 10);
		h_Pt_e[i]          = new TH1F("hPt_e_"+data_type[i], Form("Pt_e_%i_%i;P^{T}_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 2.5);
		h_theta_e[i]       = new TH1F("hTheta_e_"+data_type[i], Form("theta_e_%i_%i;#theta_{e} [rad];Counts [a.u.]", j, k), 100, 0, 40);
		h_Vz_e[i]          = new TH1F("hVz_e_"+data_type[i], Form("Vz_e_%i_%i;V_{z}^{e} [cm];Counts [a.u.]", j, k), 50, -10, 10);
		h_p_e[i]           = new TH1F("hP_e_"+data_type[i], Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 11);
		h_phi_e[i]         = new TH1F("hPhi_e_"+data_type[i], Form("phi_e_%i_%i;#phi_{e} [rad];Counts [a.u.]", j, k), 360, -180, 180);
	
		h_Z[i]             = new TH1F("hZ_"+data_type[i], Form("Z_%i_%i;Z;Counts [a.u.]", j, k), 250, 0, 1);
		h_p_pi[i]          = new TH1F("hP_pi_"+data_type[i], Form("p_pi_%i_%i;p_{#pi} [GeV];Counts [a.u.]", j, k), 250, 0, 6);
		h_Vz_pi[i]         = new TH1F("hVz_pi_"+data_type[i], Form("Vx_pi_%i_%i;V_{z}^{#pi} [cm];Counts [a.u.]", j, k), 50, -10, 10);
		h_delta_phi_pi[i]  = new TH1F("hDeltaPhi_"+data_type[i], Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 50, -360, 360);
		h_theta_pi[i]      = new TH1F("hTheta_pi_"+data_type[i], Form("theta_pi_%i_%i;#theta_{#pi} [rad];Counts [a.u.]", j, k), 100, 0, 40);
		h_phi_pi[i]        = new TH1F("hPhi_pi_"+data_type[i], Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 100, -180, 180);
		h_Pt_pi[i]         = new TH1F("hPt_pi_"+data_type[i], Form("Pt_pi_%i_%i;P^{T}_{#pi} [GeV];Counts [a.u.]", j, k), 50, 0, 1.3);
		h_Mx[i]            = new TH1F("hMx_"+data_type[i], Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 250, 0, 4);

		hQ2_W[i]		= new TH2F("hQ2_W_"+data_type[i], "; Q^{2} [GeV^{2}];W [GeV]", 250, 0, 10, 250, 1.5 ,4.5 );
		hQ2_omega[i]		= new TH2F("hQ2_omega_"+data_type[i], ";Q^{2} [GeV^{2}];#omega [GeV]", 250, 0, 10, 250, 1.5 ,10 );
		hZ_Mx[i] 		= new TH2F("hZ_Mx_"+data_type[i], ";z;M_{X} [GeV]", 50, .3, .8, 50, 1.5 ,3.5 );
		hQ2_Mx[i]		= new TH2F("hQ2_Mx_"+data_type[i], ";Q^{2} [GeV^{2}];M_{X} [GeV]", 50, 2, 8, 50, 1.5 ,3.5 );
		hXb_Mx[i]		= new TH2F("hXb_Mx_"+data_type[i], ";x_{B};M_{X} [GeV]", 50, .1, .6, 50, 1.5 ,3.5 );
		hPt_Mx[i]		= new TH2F("hPt_Mx_"+data_type[i], ";p_{#pi}^{#perp} [GeV];M_{X} [GeV]", 50, 0, 1.3, 50, 1.5 ,3.5 );
	}

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);

	//Load input tree
        TTreeReader reader_rec( chain );
	TTreeReaderValue<electron> e(reader_rec, "e");
        TTreeReaderArray<pion> pi(reader_rec, "pi");

       	int counts[2][10] = {0};

	cout<<"BEGIN EVENT LOOP\n";
	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}


		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			int idx = (int)( pi[i].getCharge() < 0 );
	


			hQ2_W[idx]->Fill( e->getQ2(), sqrt( e->getW2() ) );
			hQ2_omega[idx]->Fill( e->getQ2(), e->getOmega()  );
			
			counts[idx][0]++;
			h_W[idx]->Fill( sqrt( e->getW2() ) );
			if( sqrt(e->getW2()) < 2.5  ){continue;}
			
			counts[idx][1]++;
			h_Q2[idx]->Fill( e->getQ2() );
			if( e->getQ2() < 2 ){ continue; }

			counts[idx][2]++;
			h_y[idx]->Fill( e->getY() );
			if( e->getY() > 0.75 ){continue;}

			counts[idx][3]++;
			h_Mx[idx]->Fill( pi[i].getMx() );
			if( pi[i].getMx() < 1.7 ){ continue; }

			counts[idx][4]++;
			h_p_pi[idx]->Fill( pi[i].get3Momentum().Mag() );
			if( pi[i].get3Momentum().Mag() < 1.25 || pi[i].get3Momentum().Mag() > 5 ){continue;}

			counts[idx][5]++;
			h_theta_e[idx]->Fill( e->get3Momentum().Theta()*rad_to_deg );
			if( e->get3Momentum().Theta()*rad_to_deg < 5 || e->get3Momentum().Theta()*rad_to_deg > 40 ){continue;}
				
			counts[idx][6]++;
			h_theta_pi[idx]->Fill( pi[i].get3Momentum().Theta()*rad_to_deg );
			if( pi[i].get3Momentum().Theta()*rad_to_deg < 5 || pi[i].get3Momentum().Theta()*rad_to_deg > 40 ){continue;}
			
			counts[idx][7]++;
			h_Z[idx]->Fill( pi[i].getZ() );
			counts[idx][8]++;
		}
	}
	
	cout<<"Completed good event list... \n";

        cout<<"Writing to file\n";
        
	TFile * outFile = new TFile( out_name+".root", "RECREATE");
	outFile->cd();
	for( int i = 0; i < 2; i++ ){
		h_W[i]->Write();             
		h_Q2[i]->Write();            
		h_y[i]->Write();             
		h_theta_e[i]->Write();       
	
		h_Z[i]->Write();             
		h_p_pi[i]->Write();          
		h_theta_pi[i]->Write();      
		h_Mx[i]->Write();            

		hQ2_W[i]->Write();		
		hQ2_omega[i]->Write();		
	}

        outFile->Close();
		
	std::ofstream txtFile;
        txtFile.open(out_name + ".txt");
        txtFile<< "\t(e,e'pi+)\t#(e,e'pi-)\n";
        txtFile<< "All events\t"<<counts[0][0]<<"\t"<<counts[1][0]<<std::endl;
        txtFile<< "W2\t"<<counts[0][1]<<"\t"<<counts[1][1]<<std::endl;
        txtFile<< "Q2\t"<<counts[0][2]<<"\t"<<counts[1][2]<<std::endl;
        txtFile<< "y\t"<<counts[0][3]<<"\t"<<counts[1][3]<<std::endl;
        txtFile<< "Mx\t"<<counts[0][4]<<"\t"<<counts[1][4]<<std::endl;
        txtFile<< "p_pi\t"<<counts[0][5]<<"\t"<<counts[1][5]<<std::endl;
        txtFile<< "theta_e\t"<<counts[0][6]<<"\t"<<counts[1][6]<<std::endl;
        txtFile<< "theta_pi\t"<<counts[0][7]<<"\t"<<counts[1][7]<<std::endl;
		txtFile<< "z\t"<<counts[0][8]<<"\t"<<counts[1][8]<<std::endl;
        txtFile.close();
	

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout<<"Done. Elapsed Time : "<<elapsed.count()<<std::endl;

	return 0;
}
