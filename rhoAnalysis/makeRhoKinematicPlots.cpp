#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"
#include "analyzer.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST


using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

using namespace cutVals;
using namespace constants;

double compXf(TLorentzVector gamma, TLorentzVector pi4, double W);
TLorentzVector compPi_q(TVector3 pe, TLorentzVector q, pion pi);

enum eBkgCat { kAll = 0, kSignal = 1, kBackground = 2 };
static const int nBkgCat = 3;

static const std::string bkgCatName[nBkgCat] = {
    "all",
    "signal",
    "background"
};


int main( int argc, char** argv){

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File] [Use RICH?]  [p cut?] [Input File]\n";
		return -1;
	}

	TString out_name = argv[1];
	
	TChain * file_rec = new TChain("ePi");
	


	TString rho_norm_name = argv[2];

	cerr << "Files used: \n";
	for( int i = 3; i < argc; i++ ){	
		file_rec->Add((TString) argv[i]);
		cout<<argv[i]<<std::endl;
	}

    TFile * outFile = new TFile(out_name, "RECREATE");
	

	// Declare histograms
	TFile * rho_norms = new TFile((TString) _DATA + (TString)"/correctionFiles/" + rho_norm_name);
	TH3F * hNorms[2];
	hNorms[0] = (TH3F*)rho_norms->Get("hNorm_pip");
	hNorms[1] = (TH3F*)rho_norms->Get("hNorm_pim");
	TVector3 * bounds = (TVector3*)rho_norms->Get("bounds");
	double Mx_min = bounds->X();
	double Mx_max = bounds->Y();

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");

	cout<<"Creating Histograms\n";
	
	TH1F * h_Z[3][2][bins_Q2+1][bins_xB+1];


	TH1F * h_Mx[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_W[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_Pt_pi[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_Xb[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_Q2[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_eta[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_y[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_Vz_e[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_Vz_pi[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_p_e[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_omega[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_delta_theta[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_delta_phi[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_phi_pi[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_phi_e[3][2][bins_Q2+1][bins_xB+1];
	
	TH1F * h_theta_e[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_theta_pi[3][2][bins_Q2+1][bins_xB+1];
	TH1F* h_p_pi[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_Eta[3][2][bins_Q2+1][bins_xB+1];
	TH1F * h_Xf[3][2][bins_Q2+1][bins_xB+1];

	TH1F * h_Phi_q[3][2][bins_Q2+1][bins_xB+1];

	TH2F * hQ2_Xb[3][2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_omega[3][2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_W[3][2][bins_Q2+1][bins_xB+1];
	TH2F * hQ2_Z[3][2][bins_Q2+1][bins_xB+1];
	
	TH2F * hBeta_p[3][2][bins_Q2+1][bins_xB+1];
	TH2F * hTheta_p[3][2][bins_Q2+1][bins_xB+1];

	TString data_type[2] = {"pip", "pim"};

	for (int b = 0; b < nBkgCat; b++) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j <= bins_Q2; j++) {
				for (int k = 0; k <= bins_xB; k++) {
	
					TString tag = "_" + bkgCatName[b] + "_" + data_type[i]
								  + Form("_%i_%i", j, k);
	
					// -------------------- Electron kinematics --------------------
					h_W[b][i][j][k] =
						new TH1F("hW"+tag,
								 Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k),
								 100, 2.5, 3.7);
	
					h_Xb[b][i][j][k] =
						new TH1F("hXb"+tag,
								 Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k),
								 100, 0.15, 0.6);
	
					h_Q2[b][i][j][k] =
						new TH1F("hQ2"+tag,
								 Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k),
								 100, 0, 8);
	
					h_y[b][i][j][k] =
						new TH1F("hY"+tag,
								 Form("y_%i_%i;y;Counts [a.u.]", j, k),
								 100, 0, 1);
	
					h_omega[b][i][j][k] =
						new TH1F("hOmega"+tag,
								 Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k),
								 100, 0, 10);
	
					h_theta_e[b][i][j][k] =
						new TH1F("hTheta_e"+tag,
								 Form("theta_e_%i_%i;#theta_{e} [deg];Counts [a.u.]", j, k),
								 100, 0, 45);
	
					h_Vz_e[b][i][j][k] =
						new TH1F("hVz_e"+tag,
								 Form("Vz_e_%i_%i;V_{z}^{e} [cm];Counts [a.u.]", j, k),
								 100, -10, 10);
	
					h_p_e[b][i][j][k] =
						new TH1F("hP_e"+tag,
								 Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k),
								 100, 3, 7);
	
					h_phi_e[b][i][j][k] =
						new TH1F("hPhi_e"+tag,
								 Form("phi_e_%i_%i;#phi_{e} [deg];Counts [a.u.]", j, k),
								 360, -180, 180);
	
					// -------------------- Hadron kinematics --------------------
					
					h_delta_theta[b][i][j][k] =
						new TH1F("hDelta_theta"+tag,
								 Form("delta_theta_%i_%i;delta theta;Counts [a.u.]", j, k),
								 200, 0, 45);

					h_delta_phi[b][i][j][k] =
						new TH1F("hDelta_phi"+tag,
							Form("delta_phi_%i_%i;delta theta;Counts [a.u.]", j, k),
							720, -360, 360);			 
					
					h_Z[b][i][j][k] =
						new TH1F("hZ"+tag,
								 Form("Z_%i_%i;Z;Counts [a.u.]", j, k),
								 100, 0.3, 0.9);
	
					h_p_pi[b][i][j][k] =
						new TH1F("hP_pi"+tag,
								 Form("p_pi_%i_%i;p_{#pi} [GeV];Counts [a.u.]", j, k),
								 100, 0, 5);
	
					h_Vz_pi[b][i][j][k] =
						new TH1F("hVz_pi"+tag,
								 Form("Vz_pi_%i_%i;V_{z}^{#pi} [cm];Counts [a.u.]", j, k),
								 100, -10, 10);
	
					h_theta_pi[b][i][j][k] =
						new TH1F("hTheta_pi"+tag,
								 Form("theta_pi_%i_%i;#theta_{#pi} [deg];Counts [a.u.]", j, k),
								 100, 0, 45);
	
					h_phi_pi[b][i][j][k] =
						new TH1F("hPhi_pi"+tag,
								 Form("phi_pi_%i_%i;#phi_{#pi} [deg];Counts [a.u.]", j, k),
								 100, 0, 360);
	
					h_Pt_pi[b][i][j][k] =
						new TH1F("hPt_pi"+tag,
								 Form("Pt_pi_%i_%i;P^{T}_{#pi} [GeV];Counts [a.u.]", j, k),
								 100, 0, 1.3);
	
					h_Mx[b][i][j][k] =
						new TH1F("hMx"+tag,
								 Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k),
								 100, 1, 5);
	
					h_Phi_q[b][i][j][k] =
						new TH1F("hPhi_q"+tag,
								 Form("Phi_q_%i_%i;#phi_{q} [deg];Counts [a.u.]", j, k),
								 100, 0, 360);
	
					h_Eta[b][i][j][k] =
						new TH1F("hEta"+tag,
								 Form("Eta_%i_%i;#eta;Counts [a.u.]", j, k),
								 100, -2.5, 2.5);
	
					h_Xf[b][i][j][k] =
						new TH1F("hXf"+tag,
								 Form("Xf_%i_%i;x_{F};Counts [a.u.]", j, k),
								 100, -1, 1);
	
					// -------------------- 2D histograms --------------------
					hQ2_omega[b][i][j][k] =
						new TH2F("hQ2_omega"+tag, "",
								 100, 2, 8,
								 100, 2, 5);
	
					hQ2_Xb[b][i][j][k] =
						new TH2F("hQ2_Xb"+tag, "",
								 100, 2, 8,
								 100, 0, 0.7);
	
					hQ2_W[b][i][j][k] =
						new TH2F("hQ2_W"+tag, "",
								 100, 2, 8,
								 100, 1.5, 3.7);
	
					hQ2_Z[b][i][j][k] =
						new TH2F("hQ2_Z"+tag, "",
								 100, 2, 8,
								 100, 0.3, 1.0);
	
					hBeta_p[b][i][j][k] =
						new TH2F("hBeta_p"+tag, "",
								 100, 0, 5,
								 100, 0.99, 1.01);
	
					hTheta_p[b][i][j][k] =
						new TH2F("hTheta_p"+tag, "",
								 100, 0, 5,
								 100, 0, 40);
				}
			}
		}
	}
	

	cout<<"Beginning Event Loop\n";

	//TFile * file_rec = new TFile(in_name);
	TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion_vec(reader_rec, "isGoodPion");
	TTreeReaderArray<double> rhoWeight(reader_rec, "rhoWeight");
	TTreeReaderArray<double> rhoWeight_sym(reader_rec, "rhoWeight_sym");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");

	
	


	int event_count = 0;
	while (reader_rec.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<endl;}
                event_count++;

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



				int chargeIdx = 0;
				double radWeight = 1.;
				
				for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
					if( !isGoodPion_vec[i] || !(abs(pi[i].getPID()) == 211) ) {continue;}

					if( i == 0 ){
						//fill delta theta
					}

					bool matching = anal.applyAcceptanceMatching(pi[i], 2);
						//matching = isGoodPion[i]; }
					
					
	
					if( !matching ){ continue; }
					double M_x = pi[i].getMx();
					double pT_pi = pi[i].getPi_q().Pt();
					double p_pi = pi[i].get3Momentum().Mag();
					double theta_pi = pi[i].get3Momentum().Theta()*rad_to_deg;
					double phi_pi = pi[i].get3Momentum().Phi()*rad_to_deg;
					double Z = pi[i].getZ();
					double Vz_pi = pi[i].getVt().z() - Vz_e;
					TLorentzVector pi_q= compPi_q(e->get3Momentum(), e->getQ(), pi[i]);
					double eta = pi[i].getEta();
					double phi_q = pi_q.Phi()*rad_to_deg;
					if( phi_q < 0) phi_q += 360;
					if( phi_pi < 0) phi_pi += 360;

					double xF = compXf(e->getQ(), pi_q, W );


					int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
					int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
					int this_bin_Z = (int)( ( (Z - .3)/(1.-.3) )*bins_Z) + 1;

					int nDat = 0;
					if(*Mx_2pi > Mx_min && *Mx_2pi < Mx_max) nDat = 2;
					else if( *Mx_2pi < Mx_min) nDat = 1;
					else{continue;}
					
					for( int sym = 0; sym < 2; sym++ ){
						radWeight = (sym==0)?rhoWeight[i]:rhoWeight_sym[i];
						if( radWeight < 0 || radWeight > 100 ){continue;}
						radWeight = 0.5*radWeight;
						if( nDat == 2) radWeight *= hNorms[sym]->GetBinContent(this_bin_xB, this_bin_Q2, this_bin_Z);

						h_W[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(W, radWeight);
						h_Xb[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(xB, radWeight);
						h_Q2[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(Q2, radWeight);
						h_y[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(y, radWeight);
						h_omega[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(omega, radWeight);
						h_theta_e[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(theta_e, radWeight);
						h_Vz_e[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(Vz_e, radWeight);
						h_phi_e[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(phi_e, radWeight);
						
						if( i == 0 ){
							//fill delta theta
							TVector3 rho = pi[0].get3Momentum() + pi[1].get3Momentum();
							double delta_theta = abs(rho.Angle(e->getQ().Vect()))*rad_to_deg;
							h_delta_theta[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(delta_theta, radWeight);
							h_delta_theta[nDat][sym][0][0]->Fill(delta_theta, radWeight);

							TLorentzVector pim_q= compPi_q(e->get3Momentum(), e->getQ(), pi[1]);
							double phi_q_1 = pim_q.Phi()*rad_to_deg;
							if( phi_q_1 < 0) phi_q_1 += 360;

							h_delta_phi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( phi_q - phi_q_1);
							h_delta_phi[nDat][sym][0][0]->Fill( phi_q - phi_q_1);

						}
						h_Z[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(Z, radWeight);
						h_Mx[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(M_x, radWeight);
						h_Pt_pi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(pT_pi, radWeight);
						h_theta_pi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(theta_pi, radWeight);
						h_phi_pi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(phi_pi, radWeight);
						h_Vz_pi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(Vz_pi, radWeight);
						h_p_pi[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(p_pi, radWeight);

						h_Phi_q[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(phi_q, radWeight);
						h_Eta[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(eta, radWeight);
						h_Xf[nDat][sym][this_bin_Q2][this_bin_xB]->Fill(xF, radWeight);

						hQ2_Xb[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( Q2, xB);
						hQ2_omega[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( Q2, omega);
						hQ2_Z[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( Q2, Z);
						hQ2_W[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( xB, W);
					
						hTheta_p[nDat][sym][this_bin_Q2][this_bin_xB]->Fill( p_pi, theta_pi);

						h_W[nDat][sym][0][0]->Fill(W, radWeight);
						h_Xb[nDat][sym][0][0]->Fill(xB, radWeight);
						h_Q2[nDat][sym][0][0]->Fill(Q2, radWeight);
						h_y[nDat][sym][0][0]->Fill(y, radWeight);
						h_omega[nDat][sym][0][0]->Fill(omega, radWeight);
						h_theta_e[nDat][sym][0][0]->Fill(theta_e, radWeight);
						h_Vz_e[nDat][sym][0][0]->Fill(Vz_e, radWeight);
						h_phi_e[nDat][sym][0][0]->Fill(phi_e, radWeight);
						
						h_Z[nDat][sym][0][0]->Fill(Z, radWeight);
						h_Mx[nDat][sym][0][0]->Fill(M_x, radWeight);
						h_Pt_pi[nDat][sym][0][0]->Fill(pT_pi, radWeight);
						h_theta_pi[nDat][sym][0][0]->Fill(theta_pi, radWeight);
						h_phi_pi[nDat][sym][0][0]->Fill(phi_pi, radWeight);
						h_Vz_pi[nDat][sym][0][0]->Fill(Vz_pi, radWeight);
						h_p_pi[nDat][sym][0][0]->Fill(p_pi, radWeight);

						h_Phi_q[nDat][sym][0][0]->Fill(phi_q, radWeight);
						h_Eta[nDat][sym][0][0]->Fill(eta, radWeight);
						h_Xf[nDat][sym][0][0]->Fill(xF, radWeight);

						hQ2_Xb[nDat][sym][0][0]->Fill( Q2, xB);
						hQ2_omega[nDat][sym][0][0]->Fill( Q2, omega);
						hQ2_Z[nDat][sym][0][0]->Fill( Q2, Z);
						hQ2_W[nDat][sym][0][0]->Fill( xB, W);
					
						hTheta_p[nDat][sym][0][0]->Fill( p_pi, theta_pi);
					}
						
				}
			}
	
	//file_rec->Close();

	outFile->cd();
	auto writeSubtracted = [&](TH1* hSig, TH1* hBkg, const char* name) {
		TH1* hSub = (TH1*)hSig->Clone(name);
		hSub->Add(hBkg, -1);
		hSub->Write();
		delete hSub;
	};
	
	auto writeSubtracted2D = [&](TH2* hSig, TH2* hBkg, const char* name) {
		TH2* hSub = (TH2*)hSig->Clone(name);
		hSub->Add(hBkg, -1);
		hSub->Write();
		delete hSub;
	};

	for( int j = 0; j <= bins_Q2; j++ ){
		for( int k = 0; k <= bins_xB; k++ ){ 
			for( int i = 0; i < 2; i++ ){
				for( int nDat = 1; nDat <= 2; nDat++ ){
					h_W[nDat][i][j][k]->Write();
					h_Xb[nDat][i][j][k]->Write();
					h_Xf[nDat][i][j][k]->Write();

					h_Q2[nDat][i][j][k]->Write();

					h_Eta[nDat][i][j][k]->Write();
					h_y[nDat][i][j][k] ->Write();
					h_omega[nDat][i][j][k]->Write();
					h_theta_e[nDat][i][j][k]->Write();
					h_Vz_e[nDat][i][j][k]  ->Write();
					h_p_e[nDat][i][j][k]  ->Write();
					h_phi_e[nDat][i][j][k]->Write();
				
					h_Z[nDat][i][j][k] ->Write();
					h_p_pi[nDat][i][j][k]->Write();
					h_Vz_pi[nDat][i][j][k]->Write();
					h_theta_pi[nDat][i][j][k] ->Write();
					h_phi_pi[nDat][i][j][k]   ->Write();
					h_Phi_q[nDat][i][j][k]   ->Write();
					h_Pt_pi[nDat][i][j][k]   ->Write();
					h_Mx[nDat][i][j][k]     ->Write();
					
					hQ2_Xb[nDat][i][j][k]->Write();
					hQ2_omega[nDat][i][j][k]->Write();
					hQ2_Z[nDat][i][j][k]->Write();
					hQ2_W[nDat][i][j][k]->Write();
				
					hBeta_p[nDat][i][j][k]->Write();
					hTheta_p[nDat][i][j][k]->Write();
					
					h_delta_phi[nDat][i][j][k]->Write();
					h_delta_theta[nDat][i][j][k]->Write();

				}
				
							
				// -------------------- 1D subtracted --------------------
				writeSubtracted(h_W[1][i][j][k],      h_W[2][i][j][k],      Form("hW_sub_%i_%i_%i", i,j,k) );
				writeSubtracted(h_Xb[1][i][j][k],     h_Xb[2][i][j][k],     Form("hXb_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Xf[1][i][j][k],     h_Xf[2][i][j][k],     Form("hXf_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Q2[1][i][j][k],     h_Q2[2][i][j][k],     Form("hQ2_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Eta[1][i][j][k],    h_Eta[2][i][j][k],    Form("hEta_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_y[1][i][j][k],      h_y[2][i][j][k],      Form("hY_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_omega[1][i][j][k],  h_omega[2][i][j][k],  Form("hOmega_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_theta_e[1][i][j][k],h_theta_e[2][i][j][k],Form("hTheta_e_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Vz_e[1][i][j][k],   h_Vz_e[2][i][j][k],   Form("hVz_e_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_p_e[1][i][j][k],    h_p_e[2][i][j][k],    Form("hP_e_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_phi_e[1][i][j][k],  h_phi_e[2][i][j][k],  Form("hPhi_e_sub_%i_%i_%i", i,j,k));
		
				writeSubtracted(h_Z[1][i][j][k],      h_Z[2][i][j][k],      Form("hZ_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_p_pi[1][i][j][k],   h_p_pi[2][i][j][k],   Form("hP_pi_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Vz_pi[1][i][j][k],  h_Vz_pi[2][i][j][k],  Form("hVz_pi_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_theta_pi[1][i][j][k],h_theta_pi[2][i][j][k],Form("hTheta_pi_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_phi_pi[1][i][j][k], h_phi_pi[2][i][j][k], Form("hPhi_pi_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Phi_q[1][i][j][k],  h_Phi_q[2][i][j][k],  Form("hPhi_q_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Pt_pi[1][i][j][k],  h_Pt_pi[2][i][j][k],  Form("hPt_pi_sub_%i_%i_%i", i,j,k));
				writeSubtracted(h_Mx[1][i][j][k],     h_Mx[2][i][j][k],     Form("hMx_sub_%i_%i_%i", i,j,k));
				
				// -------------------- 2D subtracted --------------------
				//writeSubtracted2D(hQ2_Xb[1][i][j][k],     hQ2_Xb[2][i][j][k],     Form("hQ2_Xb_sub_%i_%i_%i", i,j,k));
				//writeSubtracted2D(hQ2_omega[1][i][j][k], hQ2_omega[2][i][j][k], Form("hQ2_omega_sub_%i_%i_%i", i,j,k));
				//writeSubtracted2D(hQ2_Z[1][i][j][k],      hQ2_Z[2][i][j][k],      Form("hQ2_Z_sub_%i_%i_%i", i,j,k));
				//writeSubtracted2D(hQ2_W[1][i][j][k],      hQ2_W[2][i][j][k],      Form("hQ2_W_sub_%i_%i_%i", i,j,k));
		
				//writeSubtracted2D(hTheta_p[1][i][j][k],  hTheta_p[2][i][j][k],  Form("hTheta_p_sub_%i_%i_%i", i,j,k));
					
				
				
			}
		}
	}
	outFile->Close();
}

double compXf(TLorentzVector gamma, TLorentzVector pi4, double W){
	TLorentzVector q(0, 0, gamma.Mag(), gamma.E());
	TLorentzVector targ(0,0,0, 0.938); //nucleon target
	TLorentzVector CoM_vec = (targ + gamma); //boost vector
	TVector3 CoMBoost = -1*CoM_vec.BoostVector();
	
	q = gamma;
	q.RotateZ(-gamma.Phi());
	q.RotateY(-gamma.Theta());

	q.Boost(CoMBoost);
	pi4.Boost(CoMBoost);

	TVector3 pi3 = pi4.Vect();
	TVector3 q3 = q.Vect();

	double xF = 2* pi3.Dot(q3) / (W* q3.Mag());
	return xF;
}

TLorentzVector compPi_q(TVector3 pe, TLorentzVector q, pion pi){
	pe.RotateZ(-q.Phi());
	pe.RotateY(-q.Theta());

	TLorentzVector pi_q;
	TVector3 pi3_temp = pi.get3Momentum();

	pi3_temp.RotateZ( -q.Phi()  );
    pi3_temp.RotateY( -q.Theta() );
	pi3_temp.RotateZ( -pe.Phi() );
	
	pi_q.SetVectM( pi3_temp, pi.get4Momentum().M() ); 

	return pi_q;
}

// -------------------- background subtraction --------------------

