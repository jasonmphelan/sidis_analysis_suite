#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
//#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
//#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
#include "reader.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"



using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [input k file] [input r file] [Output File]\n";
		cerr << "[Acceptance Matching Type (2,3 etc)] \n";
		cerr << "[Apply Corrections? (1 - MC, 2 - MC + pi2k, 3 - MC + pi2k + k2pi)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
	TString out_name = argv[4];

        TFile * outFile = new TFile(out_name, "RECREATE");
	//TFile * radWeightFile = new TFile("/work/clas12/users/tkutz/radgen/build/RC_graph_radgen_deuterium.root");
	//TGraph2D * rad_gen = (TGraph2D *)radWeightFile->Get("rcgr");
        
	//TFile * weightFile = new TFile( "../corrections/corrections/" + corrFileName );	

	//TH3F * accWeight_pip = (TH3F *)weightFile->Get("hAccCorrectionP");
	//TH3F * accWeight_pim = (TH3F *)weightFile->Get("hAccCorrectionM");
	//TH3F * binWeight_pip = (TH3F *)weightFile->Get("hBinMigrationP");
	//TH3F * binWeight_pim = (TH3F *)weightFile->Get("hBinMigrationM");

	// Declare histograms

	cout<<"Creating Histograms\n";
	
	TH1F * h_Z[2][bins_Q2+1][bins_xB+1];

	TH1F * h_N[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Mx[2][bins_Q2+1][bins_xB+1];
	TH1F * h_W[2][bins_Q2+1][bins_xB+1];
	TH1F * h_Pt_pi[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Xb[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Q2[2][bins_Q2+1][bins_xB+1];

	TH1F * h_eta[2][bins_Q2+1][bins_xB+1];

	TH1F * h_y[2][bins_Q2+1][bins_xB+1];
	TH1F * h_p_e[2][bins_Q2+1][bins_xB+1];

	TH1F * h_omega[2][bins_Q2+1][bins_xB+1];

	TH1F * h_phi_pi[2][bins_Q2+1][bins_xB+1];
	TH1F * h_phi_e[2][bins_Q2+1][bins_xB+1];
	
	TH1F * h_delta_phi_pi[2][bins_Q2+1][bins_xB+1];

	TH1F * h_Pt_e[2][bins_Q2+1][bins_xB+1];

	TH1F * h_theta_e[2][bins_Q2+1][bins_xB+1];
	TH1F * h_theta_pi[2][bins_Q2+1][bins_xB+1];
	TH1F* h_p_pi[2][bins_Q2+1][bins_xB+1];


	TH2F * hZ_Mx[2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_Mx[2][bins_Q2+1][bins_xB+1];
	TH2F * hXb_Mx[2][bins_Q2+1][bins_xB+1];
	TH2F * hPt_Mx[2][bins_Q2+1][bins_xB+1];

	TString data_type[2] = {"pip", "pim"};

	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){
		//int k = 0;
		/*
		h_W[0][j]             = new TH1F(Form("hW_%i", j), Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k), 50, 0, 4.2);
		h_Xb[0][j]            = new TH1F(Form("hXb_%i", j), Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k), 50, 0, 1);
		h_Q2[0][j]            = new TH1F(Form("hQ2_%i", j), Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 50, 0, 10);
		h_eta[0][j]           = new TH1F(Form("hEta_%i", j), Form("eta_%i_%i;#eta;Counts [a.u.]", j, k), 50, 0, 4.6);
		h_y[0][j]             = new TH1F(Form("hY_%i", j), Form("y_%i_%i;y;Counts [a.u.]", j, k), 50, 0, 1);
		h_omega[0][j]         = new TH1F(Form("hOmega_%i", j), Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k), 50, 0, 10);
		h_Pt_e[0][j]          = new TH1F(Form("hPt_e_%i", j), Form("Pt_e_%i_%i;P^{T}_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 2.5);
		h_theta_e[0][j]       = new TH1F(Form("hTheta_e_%i",j), Form("theta_e_%i_%i;#theta_{e} [rad];Counts [a.u.]", j, k), 50, 0, .8);
		h_Vz_e[0][j]          = new TH1F(Form("hVz_e_%i", j), Form("Vz_e_%i_%i;V_{z}^{e} [cm];Counts [a.u.]", j, k), 50, -10, 10);
		h_p_e[0][j]           = new TH1F(Form("hP_e_%i", j), Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 11);
		h_phi_e[0][j]         = new TH1F(Form("hPhi_e_%i", j), Form("phi_e_%i_%i;#phi_{e} [rad];Counts [a.u.]", j, k), 360, -180, 180);
		*/
			for( int i = 0; i < 2; i++ ){
				h_W[i][j][k]             = new TH1F("hW_"+data_type[i]+Form("_%i_%i", j, k), Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k), 50, 0, 4.2);
				h_Xb[i][j][k]            = new TH1F("hXb_"+data_type[i]+Form("_%i_%i", j, k), Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k), 50, 0, 1);
				h_Q2[i][j][k]            = new TH1F("hQ2_"+data_type[i]+Form("_%i_%i", j, k), Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 50, 0, 10);
				h_eta[i][j][k]           = new TH1F("hEta_"+data_type[i]+Form("_%i_%i", j, k), Form("eta_%i_%i;#eta;Counts [a.u.]", j, k), 50, 0, 4.6);
				h_y[i][j][k]             = new TH1F("hY_"+data_type[i]+Form("_%i_%i", j, k), Form("y_%i_%i;y;Counts [a.u.]", j, k), 50, 0, 1);
				h_omega[i][j][k]         = new TH1F("hOmega_"+data_type[i]+Form("_%i_%i", j, k), Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k), 50, 0, 10);
				h_Pt_e[i][j][k]          = new TH1F("hPt_e_"+data_type[i]+Form("_%i_%i", j, k), Form("Pt_e_%i_%i;P^{T}_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 2.5);
				h_theta_e[i][j][k]       = new TH1F("hTheta_e_"+data_type[i]+Form("_%i_%i",j, k), Form("theta_e_%i_%i;#theta_{e} [rad];Counts [a.u.]", j, k), 50, 0, 45);
				h_p_e[i][j][k]           = new TH1F("hP_e_"+data_type[i]+Form("_%i_%i", j, k), Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k), 50, 0, 11);
				h_phi_e[i][j][k]         = new TH1F("hPhi_e_"+data_type[i]+Form("_%i_%i", j, k), Form("phi_e_%i_%i;#phi_{e} [rad];Counts [a.u.]", j, k), 360, -180, 180);
			
				h_Z[i][j][k]             = new TH1F("hZ_"+data_type[i]+Form("_%i_%i", j, k), Form("Z_%i_%i;Z;Counts [a.u.]", j, k), 50, .3, .8);
				h_p_pi[i][j][k]          = new TH1F("hP_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("p_pi_%i_%i;p_{#pi} [GeV];Counts [a.u.]", j, k), 50, 0, 6);
				h_delta_phi_pi[i][j][k]  = new TH1F("hDeltaPhi_"+data_type[i]+Form("_%i_%i", j, k), Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 50, -360, 360);
				h_theta_pi[i][j][k]      = new TH1F("hTheta_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("theta_pi_%i_%i;#theta_{#pi} [rad];Counts [a.u.]", j, k), 100, 0, 45);
				h_phi_pi[i][j][k]        = new TH1F("hPhi_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 100, -180, 180);
				h_Pt_pi[i][j][k]         = new TH1F("hPt_pi_"+data_type[i]+Form("_%i_%i", j, k), Form("Pt_pi_%i_%i;P^{T}_{#pi} [GeV];Counts [a.u.]", j, k), 50, 0, 1.3);
				h_Mx[i][j][k]            = new TH1F("hMx_"+data_type[i]+Form("_%i_%i", j, k), Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 50, 0, 5);

				hZ_Mx[i][j][k]		= new TH2F("hZ_Mx_"+data_type[i]+Form("_%i_%i", j, k), "", 50, .3, .8, 50, 1.5 ,3.5 );
				hQ2_Mx[i][j][k]		= new TH2F("hQ2_Mx_"+data_type[i]+Form("_%i_%i", j, k), "", 50, 2, 8, 50, 1.5 ,3.5 );
				hXb_Mx[i][j][k]		= new TH2F("hXb_Mx_"+data_type[i]+Form("_%i_%i", j, k), "", 50, .1, .6, 50, 1.5 ,3.5 );
				hPt_Mx[i][j][k]		= new TH2F("hPt_Mx_"+data_type[i]+Form("_%i_%i", j, k), "", 50, 0, 1.3, 50, 1.5 ,3.5 );
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	TFile * file_rec = new TFile(in_name);
	TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<genElectron> e_ptr(reader_rec, "e_gen");
	TTreeReaderArray<genPion> pi_vec(reader_rec, "pi_gen");
	TTreeReaderArray<bool> isGoodPion_vec(reader_rec, "isGoodPion");


	int event_count = 0;
	while (reader_rec.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<endl;}
                event_count++;

                double Q2 = e->getQ2();
                double W = e->getW2();
                double xB = e->getXb();
                double y = e->getY();
                double omega = e->getOmega();

                double pT_e = e->get3Momentum().Pt();
                double theta_e =e->get3Momentum().Theta()*rad_to_deg 
                double p_e = e->get3Momentum().Mag();
                double phi_e = e->get3Momentum().Phi()*rad_to_pi;



		int chargeIdx = 0;
		double radWeight = 1.;
		//if(weights == 1){
		//	radWeight = rad_gen->Interpolate(xB, Q2);
		//}
		/*
                h_W[chargeIdx][0]->Fill(W, radWeight);
                h_Xb[chargeIdx][0]->Fill(xB, radWeight);
                h_Q2[chargeIdx][0]->Fill(Q2, radWeight);
                h_y[chargeIdx][0]->Fill(y, radWeight);
                h_omega[chargeIdx][0]->Fill(omega, radWeight);
                h_Pt_e[chargeIdx][0]->Fill(pT_e, radWeight);
                h_theta_e[chargeIdx][0]->Fill(theta_e*rad_to_deg, radWeight);
                h_Vz_e[chargeIdx][0]->Fill(Vz_e, radWeight);
                h_p_e[chargeIdx][0]->Fill(p_e, radWeight);
                h_phi_e[chargeIdx][0]->Fill(phi_e*rad_to_deg, radWeight);

                h_W[chargeIdx][this_bin_Q2]->Fill(W, radWeight);
                h_Xb[chargeIdx][this_bin_Q2]->Fill(xB, radWeight);
                h_Q2[chargeIdx][this_bin_Q2]->Fill(Q2, radWeight);
                h_y[chargeIdx][this_bin_Q2]->Fill(y, radWeight);
                h_omega[chargeIdx][this_bin_Q2]->Fill(omega, radWeight);
                h_Pt_e[chargeIdx][this_bin_Q2]->Fill(pT_e, radWeight);
                h_theta_e[chargeIdx][this_bin_Q2]->Fill(theta_e*rad_to_deg, radWeight);
                h_Vz_e[chargeIdx][this_bin_Q2]->Fill(Vz_e, radWeight);
                h_p_e[chargeIdx][this_bin_Q2]->Fill(p_e, radWeight);
                h_phi_e[chargeIdx][this_bin_Q2]->Fill(phi_e*rad_to_deg, radWeight);
		*/
		for( int i = 0; i < (int) ( pi_vec.end() - pi_vec.begin() ); i++ ){
			if(accMatchType == 2 && !isGoodPion_vec[i]) {continue;}
			chargeIdx = (int)(charge_vec[i] < 1);
			double M_x = pi[i].getMx();
			double pT_pi = pi[i].getPi_q().Pt();
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta*rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
			double Z = pi[i].getZ();
		
			//if( pT_pi < .5 || pT_pi > .75 ){continue;}

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
			int this_bin_Z = (int)( ( (Z - .3)/(1.-.3) )*bins_Z) + 1;
			
			
			h_W[chargeIdx][0][0]->Fill(W, radWeight);
			h_Xb[chargeIdx][0][0]->Fill(xB, radWeight);
			h_Q2[chargeIdx][0][0]->Fill(Q2, radWeight);
			h_y[chargeIdx][0][0]->Fill(y, radWeight);
			h_omega[chargeIdx][0][0]->Fill(omega, radWeight);
			h_Pt_e[chargeIdx][0][0]->Fill(pT_e, radWeight);
			h_theta_e[chargeIdx][0][0]->Fill(theta_e*rad_to_deg, radWeight);
			h_Vz_e[chargeIdx][0][0]->Fill(Vz_e, radWeight);
			h_p_e[chargeIdx][0][0]->Fill(p_e, radWeight);
			h_phi_e[chargeIdx][0][0]->Fill(phi_e*rad_to_deg, radWeight);

			h_W[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(W, radWeight);
			h_Xb[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(xB, radWeight);
			h_Q2[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Q2, radWeight);
			h_y[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(y, radWeight);
			h_omega[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(omega, radWeight);
			h_Pt_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(pT_e, radWeight);
			h_theta_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(theta_e*rad_to_deg, radWeight);
			h_p_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(p_e, radWeight);
			h_phi_e[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(phi_e*rad_to_deg, radWeight);

			h_Z[chargeIdx][0][0]->Fill(Z, radWeight);
			h_Mx[chargeIdx][0][0]->Fill(M_x, radWeight);
			h_Pt_pi[chargeIdx][0][0]->Fill(pT_pi, radWeight);
			h_theta_pi[chargeIdx][0][0]->Fill(theta_pi*rad_to_deg, radWeight);
			h_phi_pi[chargeIdx][0][0]->Fill(phi_pi*rad_to_deg, radWeight);
			h_p_pi[chargeIdx][0][0]->Fill(p_pi, radWeight);
			h_delta_phi_pi[chargeIdx][0][0]->Fill(phi_pi*rad_to_deg - phi_e*rad_to_deg, radWeight);

	
			hZ_Mx[chargeIdx][0][0]->Fill( Z, M_x, radWeight);
			hQ2_Mx[chargeIdx][0][0]->Fill( Q2, M_x, radWeight);
			hXb_Mx[chargeIdx][0][0]->Fill( xB, M_x, radWeight);
			hPt_Mx[chargeIdx][0][0]->Fill( pT_pi, M_x, radWeight);
			
			h_Z[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Z, radWeight);
			h_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(M_x, radWeight);
			h_Pt_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(pT_pi, radWeight);
			h_theta_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(theta_pi*rad_to_deg, radWeight);
			h_phi_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(phi_pi*rad_to_deg, radWeight);
			h_Vz_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(Vz_pi, radWeight);
			h_p_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(p_pi, radWeight);
			h_delta_phi_pi[chargeIdx][this_bin_Q2][this_bin_xB]->Fill(phi_pi*rad_to_deg - phi_e*rad_to_deg, radWeight);

	
			hZ_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( Z, M_x);
			hQ2_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( Q2, M_x);
			hXb_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( xB, M_x);
			hPt_Mx[chargeIdx][this_bin_Q2][this_bin_xB]->Fill( pT_pi, M_x);
			
		}
	}
	
	file_rec->Close();

	outFile->cd();
	
	for( int j = 1; j <= bins_Q2; j++ ){
		for( int k = 1; k <= bins_xB; k++ ){ 
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
				h_Pt_e[i][j][k]->Write();
				h_theta_e[i][j][k]->Write();
				h_p_e[i][j][k]  ->Write();
				h_phi_e[i][j][k]->Write();
			
				h_Z[i][j][k] ->Write();
				h_p_pi[i][j][k]->Write();
				h_delta_phi_pi[i][j][k]->Write();
				h_theta_pi[i][j][k] ->Write();
				h_phi_pi[i][j][k]   ->Write();
				h_Pt_pi[i][j][k]   ->Write();
				h_Mx[i][j][k]     ->Write();
				
				hZ_Mx[i][j][k]->Write();
				hQ2_Mx[i][j][k]->Write();
				hXb_Mx[i][j][k]->Write();
				hPt_Mx[i][j][k]->Write();

			}
		}
	}
	outFile->Close();
}
