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

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input path] [Output File] [# of input files] [Beam Energy] \n";
		return -1;
	}
	cerr << "Files used: " << argv[2] << "\nnFiles " << atoi(argv[3]) << "\n";

	TString in_name = argv[1];
    	TString out_name = argv[2];
    	int nFiles = atoi(argv[3]);
    	double EBeam = atof(argv[4]);
        
	TFile * outFile = new TFile(out_name +".root", "RECREATE");

	reader skimReader;
	skimReader.setNumFiles( nFiles);
	skimReader.setRunType( 1 );
	skimReader.setEnergy( EBeam );

	TChain * chain = new TChain("ePi");
	skimReader.getRunSkimsByName(chain, in_name);

        //TFile * file_rec = new TFile(in_name, "UPDATE");
	
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);
	
	TH1F * h_Z[2][bins_Q2+1][bins_xB+1];


	TH1F * h_Mx[2][bins_Q2+1][bins_xB+1];
	TH1F * h_W[2][bins_Q2+1][bins_xB+1];
	TH1F * h_Pt_pi[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Xb[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Q2[2][bins_Q2+1][bins_xB+1];

	TH1F * h_eta[2][bins_Q2+1][bins_xB+1];

	TH1F * h_y[2][bins_Q2+1][bins_xB+1];
	TH1F * h_Vz_e[2][bins_Q2+1][bins_xB+1];
	TH1F * h_Vz_pi[2][bins_Q2+1][bins_xB+1];
	TH1F * h_p_e[2][bins_Q2+1][bins_xB+1];

	TH1F * h_omega[2][bins_Q2+1][bins_xB+1];

	TH1F * h_phi_pi[2][bins_Q2+1][bins_xB+1];
	TH1F * h_phi_e[2][bins_Q2+1][bins_xB+1];
	
	TH1F * h_theta_e[2][bins_Q2+1][bins_xB+1];
	TH1F * h_theta_pi[2][bins_Q2+1][bins_xB+1];
	TH1F* h_p_pi[2][bins_Q2+1][bins_xB+1];


	TH2F * hQ2_omega[2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_W[2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_Z[2][bins_Q2+1][bins_xB+1];

	TString data_type[2] = {"pip", "pim"};

	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){
			for( int i = 0; i < 2; i++ ){
				h_W[i][j][k]             = new TH1F("hW_"+data_type[i]+Form("_%i_%i", j, k), Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k), 250, 2.5, 3.7);
				h_Xb[i][j][k]            = new TH1F("hXb_"+data_type[i]+Form("_%i_%i", j, k), Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k), 250, 0.15, .6);
				h_Q2[i][j][k]            = new TH1F("hQ2_"+data_type[i]+Form("_%i_%i", j, k), Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 250, 0, 8);
				h_eta[i][j][k]           = new TH1F("hEta_"+data_type[i]+Form("_%i_%i", j, k), Form("eta_%i_%i;#eta;Counts [a.u.]", j, k), 250, 0.5, 4.5);
				h_y[i][j][k]             = new TH1F("hY_"+data_type[i]+Form("_%i_%i", j, k), Form("y_%i_%i;y;Counts [a.u.]", j, k), 250, 0, 1);
				h_omega[i][j][k]         = new TH1F("hOmega_"+data_type[i]+Form("_%i_%i", j, k), Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k), 250, 0, 10);
				h_theta_e[i][j][k]       = new TH1F("hTheta_e_"+data_type[i]+Form("_%i_%i",j, k), Form("theta_e_%i_%i;#theta_{e} [rad];Counts [a.u.]", j, k), 250, 0, 45);
				h_Vz_e[i][j][k]          = new TH1F("hVz_e_"+data_type[i]+Form("_%i_%i", j, k), Form("Vz_e_%i_%i;V_{z}^{e} [cm];Counts [a.u.]", j, k), 250, -10, 10);
				h_p_e[i][j][k]           = new TH1F("hP_e_"+data_type[i]+Form("_%i_%i", j, k), Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k), 250, 3, 7);
				h_phi_e[i][j][k]         = new TH1F("hPhi_e_"+data_type[i]+Form("_%i_%i", j, k), Form("phi_e_%i_%i;#phi_{e} [rad];Counts [a.u.]", j, k), 360, -180, 180);
			
				h_Z[i][j][k]             = new TH1F("hZ_"+data_type[i]+Form("_%i_%i", j, k), Form("Z_%i_%i;Z;Counts [a.u.]", j, k), 250, .3, .9);
				h_p_pi[i][j][k]          = new TH1F("hP_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("p_pi_%i_%i;p_{#pi} [GeV];Counts [a.u.]", j, k), 250, 0, 5);
				h_Vz_pi[i][j][k]         = new TH1F("hVz_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("Vx_pi_%i_%i;V_{z}^{#pi} [cm];Counts [a.u.]", j, k), 250, -10, 10);
				h_theta_pi[i][j][k]      = new TH1F("hTheta_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("theta_pi_%i_%i;#theta_{#pi} [rad];Counts [a.u.]", j, k), 250, 0, 45);
				h_phi_pi[i][j][k]        = new TH1F("hPhi_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 250, -180, 180);
				h_Pt_pi[i][j][k]         = new TH1F("hPt_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("Pt_pi_%i_%i;P^{T}_{#pi} [GeV];Counts [a.u.]", j, k), 250, 0, 1.3);
				h_Mx[i][j][k]            = new TH1F("hMx_"+data_type[i]+Form("_%i_%i", j, k), Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 250, 1, 5);

				hQ2_omega[i][j][k]		= new TH2F("hQ2_omega_"+data_type[i]+Form("_%i_%i", j, k), "", 250, 2, 8, 250, 2 ,5 );
				hQ2_W[i][j][k]		= new TH2F("hQ2_W_"+data_type[i]+Form("_%i_%i", j, k), "", 250, 2, 8, 250, 1.5 ,3.7 );
				hQ2_Z[i][j][k]		= new TH2F("hQ2_Z_"+data_type[i]+Form("_%i_%i", j, k), "", 250, 2, 8, 250, .3 ,1 );
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	//Load input tree
        //TTreeReader reader_rec("ePi", file_rec);
    	TTreeReader reader_rec( chain );
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	
	int event_total = reader_rec.GetEntries();
	double counts[2][20] = {0}; //[charge][cut level]

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

		int radWeight = 1;
		int chargeIdx;
                double Q2 = e->getQ2();
                double W = sqrt(e->getW2());
                double xB = e->getXb();
                double y = e->getY();
                double omega = e->getOmega();

                double pT_e = e->get3Momentum().Pt();
                double theta_e =e->get3Momentum().Theta()*rad_to_deg;
                double p_e = e->get3Momentum().Mag();
                double phi_e = e->get3Momentum().Phi()*rad_to_deg;

                double Vz_e = e->getVt().z();



		//if( !anal.applyElectronKinematicCuts( *e ) ){ continue; }

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
               
                	//if(!anal.applyPionKinematicCuts(pi[i])){ continue; }
			chargeIdx = (int)(pi[i].getCharge() < 1);
			double M_x = pi[i].getMx();
			double pT_pi = pi[i].getPi_q().Pt();
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()*rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
			double Z = pi[i].getZ();
			double Vz_pi = pi[i].getVt().z() - Vz_e;
		
			//if( pT_pi < .5 || pT_pi > .75 ){continue;}

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
			int this_bin_Z = (int)( ( (Z - .3)/(1.-.3) )*bins_Z) + 1;
			
			
			h_W[chargeIdx][0][0]->Fill(W, radWeight);
			h_Xb[chargeIdx][0][0]->Fill(xB, radWeight);
			h_Q2[chargeIdx][0][0]->Fill(Q2, radWeight);
			h_y[chargeIdx][0][0]->Fill(y, radWeight);
			h_omega[chargeIdx][0][0]->Fill(omega, radWeight);
			h_theta_e[chargeIdx][0][0]->Fill(theta_e*rad_to_deg, radWeight);
			h_Vz_e[chargeIdx][0][0]->Fill(Vz_e, radWeight);
			h_p_e[chargeIdx][0][0]->Fill(p_e, radWeight);
			h_phi_e[chargeIdx][0][0]->Fill(phi_e*rad_to_deg, radWeight);

			h_W[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(W, radWeight);
			h_Xb[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(xB, radWeight);
			h_Q2[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Q2, radWeight);
			h_y[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(y, radWeight);
			h_omega[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(omega, radWeight);
			h_theta_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(theta_e*rad_to_deg, radWeight);
			h_Vz_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Vz_e, radWeight);
			h_phi_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(phi_e*rad_to_deg, radWeight);

			h_Z[chargeIdx][0][0]->Fill(Z, radWeight);
			h_Mx[chargeIdx][0][0]->Fill(M_x, radWeight);
			h_Pt_pi[chargeIdx][0][0]->Fill(pT_pi, radWeight);
			h_theta_pi[chargeIdx][0][0]->Fill(theta_pi*rad_to_deg, radWeight);
			h_phi_pi[chargeIdx][0][0]->Fill(phi_pi*rad_to_deg, radWeight);
			h_Vz_pi[chargeIdx][0][0]->Fill(Vz_pi - Vz_e, radWeight);
			h_p_pi[chargeIdx][0][0]->Fill(p_pi, radWeight);
	
			hQ2_omega[chargeIdx][0][0]->Fill( Q2, omega, radWeight);
			hQ2_Z[chargeIdx][0][0]->Fill( Q2, Z, radWeight);
			hQ2_W[chargeIdx][0][0]->Fill( Q2, W, radWeight);
			
			h_Z[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Z, radWeight);
			h_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(M_x, radWeight);
			h_Pt_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(pT_pi, radWeight);
			h_theta_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(theta_pi*rad_to_deg, radWeight);
			h_phi_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(phi_pi*rad_to_deg, radWeight);
			h_Vz_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Vz_pi, radWeight);
			h_p_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(p_pi, radWeight);

	
			hQ2_omega[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( Q2, omega);
			hQ2_Z[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( Q2, Z);
			hQ2_W[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( xB, W);
			

		}
	
	}
	
	std::ofstream txtFile;
	txtFile.open(out_name + ".txt");
	txtFile<< "\t(e,e'pi+)\t#(e,e'pi-)\n";
	txtFile<< "All tracks\t"<<counts[0][0]<<"\t"<<counts[1][0]<<std::endl;
	txtFile<< "Event Builder\t"<<counts[0][1]<<"\t"<<counts[1][1]<<std::endl;
	txtFile<< "Electron DC Fiducials\t"<<counts[0][2]<<"\t"<<counts[1][2]<<std::endl;
	txtFile<< "PCAL WV\t"<<counts[0][3]<<"\t"<<counts[1][3]<<std::endl;
	txtFile<< "PCAL Edep\t"<<counts[0][4]<<"\t"<<counts[1][4]<<std::endl;
	txtFile<< "SF Cuts\t"<<counts[0][5]<<"\t"<<counts[1][5]<<std::endl;
	txtFile<< "SF Correlation\t"<<counts[0][6]<<"\t"<<counts[1][6]<<std::endl;
	txtFile<< "Electron Vertex\t"<<counts[0][7]<<"\t"<<counts[1][7]<<std::endl;
	txtFile<< "Pion DC Fiducials\t"<<counts[0][8]<<"\t"<<counts[1][8]<<std::endl;
	txtFile<< "Pion Vertex\t"<<counts[0][9]<<"\t"<<counts[1][9]<<std::endl;
	txtFile<< "Chi2\t"<<counts[0][10]<<"\t"<<counts[1][10]<<std::endl;
	txtFile.close();


	outFile->cd();
	
	for( int j = 0; j <= 0; j++ ){
		for( int k = 0; k <= 0; k++ ){ 
			/*
			h_W[0][j]->Write();
			h_Xb[0][j]->Write();

			h_Q2[0][j]->Write();

			h_eta[0][j]->Write();
			h_y[0][j] ->Write();
			h_omega[0][j]->Write();
			h_Pt_e[0][j]->Write();
			h_theta_e[0][j]->Write();
			h_Vz_e[0][j]  ->Write();
			h_p_e[0][j]  ->Write();
			h_phi_e[0][j]->Write();
	*/		
			for( int i = 0; i < 2; i++ ){
				h_W[i][j][k]->Write();
				h_Xb[i][j][k]->Write();

				h_Q2[i][j][k]->Write();

				h_eta[i][j][k]->Write();
				h_y[i][j][k] ->Write();
				h_omega[i][j][k]->Write();
				h_theta_e[i][j][k]->Write();
				h_Vz_e[i][j][k]  ->Write();
				h_p_e[i][j][k]  ->Write();
				h_phi_e[i][j][k]->Write();
			
				h_Z[i][j][k] ->Write();
				h_p_pi[i][j][k]->Write();
				h_Vz_pi[i][j][k]->Write();
				h_theta_pi[i][j][k] ->Write();
				h_phi_pi[i][j][k]   ->Write();
				h_Pt_pi[i][j][k]   ->Write();
				h_Mx[i][j][k]     ->Write();
				
				hQ2_omega[i][j][k]->Write();
				hQ2_Z[i][j][k]->Write();
				hQ2_W[i][j][k]->Write();

			}
		}
	}
	outFile->Close();

	return 0;
}
