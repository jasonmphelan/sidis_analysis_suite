#include <iostream>
#include <cmath>
#include <fstream>
#include <memory>
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

static const int MAX_VAR_BINS = 20;

double getVarVal(TString var, double xB, double Q2, double Z,
                 double pT_pi, double M_x, double xF,
                 double eta, double p_pi, double theta_pi, double phi_pi){
	if( var == "Z"   || var == "z"  ) return Z;
	if( var == "pT"  || var == "Pt" ) return pT_pi;
	if( var == "Mx"  || var == "M_x") return M_x;
	if( var == "xF"                 ) return xF;
	if( var == "eta"                ) return eta;
	if( var == "p_pi"               ) return p_pi;
	if( var == "theta_pi"           ) return theta_pi;
	if( var == "phi_pi"             ) return phi_pi;
	if( var == "xB"                 ) return xB;
	if( var == "Q2"                 ) return Q2;
	return 0;
}

int main( int argc, char** argv){

	if( argc < 10 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File] [isSim (0=data,1=sim)] [accCorr (0=off,1=on)] [corrFitFile] [var_name] [bins_var] [var_min] [var_max] [Input File...]\n";
		return -1;
	}

	TString out_name = argv[1];

	TChain * file_rec = new TChain("ePi");


	int acc_match       = atoi(argv[2]);   // 0 = data, 1 = simulation
	int     doAccCorr   = atoi(argv[3]);   // 0 = no correction, 1 = apply acceptance correction
	TString corrFitFile = argv[4];         // e.g. corrections_fit.root or corrections_10.4_fit.root
	TString var_name    = argv[5];         // extra binning variable, or "null"
	int     bins_var    = atoi(argv[6]);   // number of extra variable bins (0 = no extra binning)
	double  var_min     = atof(argv[7]);
	double  var_max     = atof(argv[8]);

	cout << "isSim    : " << (acc_match    == 1 ? "simulation"   : "data") << "\n";
	cout << "accCorr  : " << (doAccCorr == 1 ? "on"           : "off")  << "\n";
	cout << "corrFile : " << corrFitFile << "\n";
	cout << "var_name : " << var_name << "  bins_var=" << bins_var
	     << "  [" << var_min << ", " << var_max << "]\n";

	cerr << "Files used: \n";
	for( int i = 9; i < argc; i++ ){
		file_rec->Add((TString) argv[i]);
		cout<<argv[i]<<std::endl;
	}

    TFile * outFile = new TFile(out_name, "RECREATE");

	correctionTools corr(bins_var > 0 ? 4 : 2);
	if( doAccCorr ){
		corr.setN4dBins(bins_var);
		corr.setWeightFitName( corrFitFile );
		corr.loadContinuousCorrections();
		cout << "Loaded continuous correction fits\n";
	}

	// Declare histograms
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//isSim);//runType);
	anal.loadSamplingFractionParams();
	anal.loadMatchingFunctions("matchCut2D_10.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));


	cout<<"Creating Histograms\n";

	TH1F * h_Z[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_Mx[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_W[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_Pt_pi[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_Xb[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_Q2[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_eta[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_y[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_Vz_e[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_Vz_pi[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_p_e[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_omega[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_phi_pi[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_phi_e[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_theta_e[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_theta_pi[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F* h_p_pi[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_Eta[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH1F * h_Xf[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH1F * h_Phi_q[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH2F * hQ2_Xb[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH2F * hQ2_omega[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH2F * hQ2_W[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH2F * hQ2_Z[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TH2F * hBeta_p[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];
	TH2F * hTheta_p[2][bins_Q2+1][bins_xB+1][MAX_VAR_BINS+1];

	TString data_type[2] = {"pip", "pim"};

	for( int v = 0; v <= bins_var; v++ ){
		TString vsuf = (v > 0) ? Form("_v%d", v) : "";
		for( int j = 0; j <= bins_Q2; j++ ){
			for( int k = 0; k <= bins_xB; k++ ){
				for( int i = 0; i < 2; i++ ){
					h_W[i][j][k][v]             = new TH1F("hW_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("W_%i_%i;W [GeV];Counts [a.u.]", j, k), 100, 2.5, 3.7);
					h_Xb[i][j][k][v]            = new TH1F("hXb_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("xB_%i_%i;x_{B};Counts [a.u.]", j, k), 100, 0.15, .6);
					h_Q2[i][j][k][v]            = new TH1F("hQ2_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Q2_%i_%i;Q^{2} [GeV^{2}];Counts [a.u.]", j, k), 100, 0, 8);
					h_y[i][j][k][v]             = new TH1F("hY_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("y_%i_%i;y;Counts [a.u.]", j, k), 100, 0, 1);
					h_omega[i][j][k][v]         = new TH1F("hOmega_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("omega_%i_%i;#omega [GeV];Counts [a.u.]", j, k), 100, 0, 10);
					h_theta_e[i][j][k][v]       = new TH1F("hTheta_e_"+data_type[i]+Form("_%i_%i",j, k)+vsuf, Form("theta_e_%i_%i;#theta_{e} [rad];Counts [a.u.]", j, k), 100, 0, 45);
					h_Vz_e[i][j][k][v]          = new TH1F("hVz_e_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Vz_e_%i_%i;V_{z}^{e} [cm];Counts [a.u.]", j, k), 100, -10, 10);
					h_p_e[i][j][k][v]           = new TH1F("hP_e_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("p_e_%i_%i;p_{e} [GeV];Counts [a.u.]", j, k), 100, 3, 7);
					h_phi_e[i][j][k][v]         = new TH1F("hPhi_e_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("phi_e_%i_%i;#phi_{e} [rad];Counts [a.u.]", j, k), 360, -180, 180);

					h_Z[i][j][k][v]             = new TH1F("hZ_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Z_%i_%i;Z;Counts [a.u.]", j, k), 100, .3, .9);
					h_p_pi[i][j][k][v]          = new TH1F("hP_pi_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("p_pi_%i_%i;p_{#pi} [GeV];Counts [a.u.]", j, k), 100, 0, 5);
					h_Vz_pi[i][j][k][v]         = new TH1F("hVz_pi_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Vx_pi_%i_%i;V_{z}^{#pi} [cm];Counts [a.u.]", j, k), 100, -10, 10);
					h_theta_pi[i][j][k][v]      = new TH1F("hTheta_pi_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("theta_pi_%i_%i;#theta_{#pi} [rad];Counts [a.u.]", j, k), 100, 0, 45);
					h_phi_pi[i][j][k][v]        = new TH1F("hPhi_pi_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("phi_pi_%i_%i;#phi_{#pi} [rad];Counts [a.u.]", j, k), 100, 0, 360);
					h_Pt_pi[i][j][k][v]         = new TH1F("hPt_pi_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Pt_pi_%i_%i;P^{T}_{#pi} [GeV];Counts [a.u.]", j, k), 100, 0, 1.3);
					h_Mx[i][j][k][v]            = new TH1F("hMx_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 100, 1, 5);
					h_Phi_q[i][j][k][v]         = new TH1F("hPhi_q_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 100, 0, 360);
					h_Eta[i][j][k][v]           = new TH1F("hEta_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 100, 1, 5);
					h_Xf[i][j][k][v]            = new TH1F("hXf_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, Form("Mx_%i_%i;M_{x} [GeV];Counts [a.u.]", j, k), 100, -1, 1);

					hQ2_omega[i][j][k][v]       = new TH2F("hQ2_omega_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 2, 8, 100, 2 ,5 );
					hQ2_Xb[i][j][k][v]         = new TH2F("hQ2_Xb_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 2, 8, 100, 0 ,0.7 );
					hQ2_W[i][j][k][v]           = new TH2F("hQ2_W_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 2, 8, 100, 1.5 ,3.7 );
					hQ2_Z[i][j][k][v]           = new TH2F("hQ2_Z_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 2, 8, 100, .3 ,1 );

					hBeta_p[i][j][k][v]         = new TH2F("hBeta_p_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 0, 5, 100, .99 ,1.01 );
					hTheta_p[i][j][k][v]        = new TH2F("hTheta_p_"+data_type[i]+Form("_%i_%i", j, k)+vsuf, "", 100, 0, 10, 100, 0 , 40 );
				}
			}
		}
	}

	cout<<"Beginning Event Loop\n";

	TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	// target branch only present in simulation trees; skip reader setup for data
	std::unique_ptr<TTreeReaderValue<int>> sim_target;

	int event_count = 0;
	while (reader_rec.Next()) {
                if( event_count == 10000000) break;
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

		if( !anal.applyElectronDetectorCuts( *e ) ) continue;
		if( !anal.applyElectronKinematicCuts( *e ) ) continue;
		if( anal.applyAcceptanceMap( e->get3Momentum().Mag(), rad_to_deg*e->get3Momentum().Phi(), rad_to_deg*e->get3Momentum().Theta(), 0 ) < 0 ) continue;


		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			if( !(abs(pi[i].getPID()) == 211) ) {continue;}
			if( !anal.applyPionDetectorCuts(pi[i], *e) ) continue;
			if( !anal.applyPionKinematicCuts(pi[i]) ) continue;
			if( acc_match && !anal.applyAcceptanceMatching(pi[i], 2) ) continue;

			chargeIdx = (int)(pi[i].getCharge() < 1);

			if(anal.applyAcceptanceMap( pi[i].get3Momentum().Mag(), rad_to_deg*pi[i].get3Momentum().Phi(), rad_to_deg*pi[i].get3Momentum().Theta(), 1+chargeIdx ) < 0 ) continue;

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

			double beta = 0;

			int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*bins_Q2) + 1;
			int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*bins_xB) + 1;
			int this_bin_Z = (int)( ( (Z - .3)/(1.-.3) )*bins_Z) + 1;

			// ── extra variable bin (must be computed before correction lookup) ────
			int this_bin_var = -1;
			if( var_name != "null" && bins_var > 0 ){
				double var_val = getVarVal(var_name, xB, Q2, Z, pT_pi, M_x, xF, eta, p_pi, theta_pi, phi_pi);
				int bv = (int)( ((var_val - var_min)/(var_max - var_min)) * bins_var );
				if( bv >= 0 && bv < bins_var ) this_bin_var = bv + 1;
			}
			// set4dBin uses 0-based index; -1 if no valid bin (falls back to 3D fit)
			corr.set4dBin( this_bin_var >= 1 ? this_bin_var - 1 : -1 );

			corr.setKinematics( xB, Q2, Z, p_pi );
			double accW = doAccCorr ? corr.getCorrectionFactor( 2, chargeIdx ) : 1.0;
			if( doAccCorr && accW <= 0 ) continue;
			const double w = radWeight * accW;

			// ── inclusive [0][0] ───────────────────────────────────────────────
			hTheta_p[chargeIdx][0][0][0]->Fill( e->get3Momentum().Mag(), e->get3Momentum().Theta()*rad_to_deg, w);

			h_W[chargeIdx][0][0][0]->Fill(W, w);
			h_Xb[chargeIdx][0][0][0]->Fill(xB, w);
			h_Q2[chargeIdx][0][0][0]->Fill(Q2, w);
			h_y[chargeIdx][0][0][0]->Fill(y, w);
			h_omega[chargeIdx][0][0][0]->Fill(omega, w);
			h_theta_e[chargeIdx][0][0][0]->Fill(theta_e, w);
			h_Vz_e[chargeIdx][0][0][0]->Fill(Vz_e, w);
			h_p_e[chargeIdx][0][0][0]->Fill(p_e, w);
			h_phi_e[chargeIdx][0][0][0]->Fill(phi_e, w);
			h_Z[chargeIdx][0][0][0]->Fill(Z, w);
			h_Mx[chargeIdx][0][0][0]->Fill(M_x, w);
			h_Pt_pi[chargeIdx][0][0][0]->Fill(pT_pi, w);
			h_theta_pi[chargeIdx][0][0][0]->Fill(theta_pi, w);
			h_phi_pi[chargeIdx][0][0][0]->Fill(phi_pi, w);
			h_Vz_pi[chargeIdx][0][0][0]->Fill(Vz_pi - Vz_e, w);
			h_p_pi[chargeIdx][0][0][0]->Fill(p_pi, w);
			h_Phi_q[chargeIdx][0][0][0]->Fill(phi_q, w);
			h_Eta[chargeIdx][0][0][0]->Fill(eta, w);
			h_Xf[chargeIdx][0][0][0]->Fill(xF, w);
			hQ2_Xb[chargeIdx][0][0][0]->Fill( Q2, xB, w);
			hQ2_omega[chargeIdx][0][0][0]->Fill( Q2, omega, w);
			hQ2_Z[chargeIdx][0][0][0]->Fill( Q2, Z, w);
			hQ2_W[chargeIdx][0][0][0]->Fill( Q2, W, w);

			// ── Q2 + xB binned [q][x] ─────────────────────────────────────────
			h_W[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(W, w);
			h_Xb[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(xB, w);
			h_Q2[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(Q2, w);
			h_y[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(y, w);
			h_omega[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(omega, w);
			h_theta_e[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(theta_e, w);
			h_Vz_e[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(Vz_e, w);
			h_phi_e[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(phi_e, w);
			h_Z[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(Z, w);
			h_Mx[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(M_x, w);
			h_Pt_pi[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(pT_pi, w);
			h_theta_pi[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(theta_pi, w);
			h_phi_pi[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(phi_pi, w);
			h_Vz_pi[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(Vz_pi, w);
			h_p_pi[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(p_pi, w);
			h_Phi_q[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(phi_q, w);
			h_Eta[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(eta, w);
			h_Xf[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill(xF, w);
			hQ2_Xb[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill( Q2, xB, w);
			hQ2_omega[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill( Q2, omega, w);
			hQ2_Z[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill( Q2, Z, w);
			hQ2_W[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill( xB, W, w);
			hBeta_p[chargeIdx][this_bin_Q2][this_bin_xB][0]->Fill( p_pi, beta, w);

			// ── Q² fixed, xB inclusive [q][0] ─────────────────────────────────
			h_W[chargeIdx][this_bin_Q2][0][0]->Fill(W, w);
			h_Xb[chargeIdx][this_bin_Q2][0][0]->Fill(xB, w);
			h_Q2[chargeIdx][this_bin_Q2][0][0]->Fill(Q2, w);
			h_y[chargeIdx][this_bin_Q2][0][0]->Fill(y, w);
			h_omega[chargeIdx][this_bin_Q2][0][0]->Fill(omega, w);
			h_theta_e[chargeIdx][this_bin_Q2][0][0]->Fill(theta_e, w);
			h_Vz_e[chargeIdx][this_bin_Q2][0][0]->Fill(Vz_e, w);
			h_p_e[chargeIdx][this_bin_Q2][0][0]->Fill(p_e, w);
			h_phi_e[chargeIdx][this_bin_Q2][0][0]->Fill(phi_e, w);
			h_Z[chargeIdx][this_bin_Q2][0][0]->Fill(Z, w);
			h_Mx[chargeIdx][this_bin_Q2][0][0]->Fill(M_x, w);
			h_Pt_pi[chargeIdx][this_bin_Q2][0][0]->Fill(pT_pi, w);
			h_theta_pi[chargeIdx][this_bin_Q2][0][0]->Fill(theta_pi, w);
			h_phi_pi[chargeIdx][this_bin_Q2][0][0]->Fill(phi_pi, w);
			h_Vz_pi[chargeIdx][this_bin_Q2][0][0]->Fill(Vz_pi, w);
			h_p_pi[chargeIdx][this_bin_Q2][0][0]->Fill(p_pi, w);
			h_Phi_q[chargeIdx][this_bin_Q2][0][0]->Fill(phi_q, w);
			h_Eta[chargeIdx][this_bin_Q2][0][0]->Fill(eta, w);
			h_Xf[chargeIdx][this_bin_Q2][0][0]->Fill(xF, w);
			hQ2_Xb[chargeIdx][this_bin_Q2][0][0]->Fill( Q2, xB, w);
			hQ2_omega[chargeIdx][this_bin_Q2][0][0]->Fill( Q2, omega, w);
			hQ2_Z[chargeIdx][this_bin_Q2][0][0]->Fill( Q2, Z, w);
			hQ2_W[chargeIdx][this_bin_Q2][0][0]->Fill( xB, W, w);

			// ── xB fixed, Q² inclusive [0][x] ─────────────────────────────────
			h_W[chargeIdx][0][this_bin_xB][0]->Fill(W, w);
			h_Xb[chargeIdx][0][this_bin_xB][0]->Fill(xB, w);
			h_Q2[chargeIdx][0][this_bin_xB][0]->Fill(Q2, w);
			h_y[chargeIdx][0][this_bin_xB][0]->Fill(y, w);
			h_omega[chargeIdx][0][this_bin_xB][0]->Fill(omega, w);
			h_theta_e[chargeIdx][0][this_bin_xB][0]->Fill(theta_e, w);
			h_Vz_e[chargeIdx][0][this_bin_xB][0]->Fill(Vz_e, w);
			h_p_e[chargeIdx][0][this_bin_xB][0]->Fill(p_e, w);
			h_phi_e[chargeIdx][0][this_bin_xB][0]->Fill(phi_e, w);
			h_Z[chargeIdx][0][this_bin_xB][0]->Fill(Z, w);
			h_Mx[chargeIdx][0][this_bin_xB][0]->Fill(M_x, w);
			h_Pt_pi[chargeIdx][0][this_bin_xB][0]->Fill(pT_pi, w);
			h_theta_pi[chargeIdx][0][this_bin_xB][0]->Fill(theta_pi, w);
			h_phi_pi[chargeIdx][0][this_bin_xB][0]->Fill(phi_pi, w);
			h_Vz_pi[chargeIdx][0][this_bin_xB][0]->Fill(Vz_pi, w);
			h_p_pi[chargeIdx][0][this_bin_xB][0]->Fill(p_pi, w);
			h_Phi_q[chargeIdx][0][this_bin_xB][0]->Fill(phi_q, w);
			h_Eta[chargeIdx][0][this_bin_xB][0]->Fill(eta, w);
			h_Xf[chargeIdx][0][this_bin_xB][0]->Fill(xF, w);
			hQ2_Xb[chargeIdx][0][this_bin_xB][0]->Fill( Q2, xB, w);
			hQ2_omega[chargeIdx][0][this_bin_xB][0]->Fill( Q2, omega, w);
			hQ2_Z[chargeIdx][0][this_bin_xB][0]->Fill( Q2, Z, w);
			hQ2_W[chargeIdx][0][this_bin_xB][0]->Fill( xB, W, w);

			// ── extra variable bin [v] ─────────────────────────────────────────
			if( this_bin_var >= 1 ){
				const int v = this_bin_var;

				hTheta_p[chargeIdx][0][0][v]->Fill( e->get3Momentum().Mag(), e->get3Momentum().Theta()*rad_to_deg, w);
				h_W[chargeIdx][0][0][v]->Fill(W, w);
				h_Xb[chargeIdx][0][0][v]->Fill(xB, w);
				h_Q2[chargeIdx][0][0][v]->Fill(Q2, w);
				h_y[chargeIdx][0][0][v]->Fill(y, w);
				h_omega[chargeIdx][0][0][v]->Fill(omega, w);
				h_theta_e[chargeIdx][0][0][v]->Fill(theta_e, w);
				h_Vz_e[chargeIdx][0][0][v]->Fill(Vz_e, w);
				h_p_e[chargeIdx][0][0][v]->Fill(p_e, w);
				h_phi_e[chargeIdx][0][0][v]->Fill(phi_e, w);
				h_Z[chargeIdx][0][0][v]->Fill(Z, w);
				h_Mx[chargeIdx][0][0][v]->Fill(M_x, w);
				h_Pt_pi[chargeIdx][0][0][v]->Fill(pT_pi, w);
				h_theta_pi[chargeIdx][0][0][v]->Fill(theta_pi, w);
				h_phi_pi[chargeIdx][0][0][v]->Fill(phi_pi, w);
				h_Vz_pi[chargeIdx][0][0][v]->Fill(Vz_pi - Vz_e, w);
				h_p_pi[chargeIdx][0][0][v]->Fill(p_pi, w);
				h_Phi_q[chargeIdx][0][0][v]->Fill(phi_q, w);
				h_Eta[chargeIdx][0][0][v]->Fill(eta, w);
				h_Xf[chargeIdx][0][0][v]->Fill(xF, w);
				hQ2_Xb[chargeIdx][0][0][v]->Fill( Q2, xB, w);
				hQ2_omega[chargeIdx][0][0][v]->Fill( Q2, omega, w);
				hQ2_Z[chargeIdx][0][0][v]->Fill( Q2, Z, w);
				hQ2_W[chargeIdx][0][0][v]->Fill( Q2, W, w);

				h_W[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(W, w);
				h_Xb[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(xB, w);
				h_Q2[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(Q2, w);
				h_y[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(y, w);
				h_omega[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(omega, w);
				h_theta_e[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(theta_e, w);
				h_Vz_e[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(Vz_e, w);
				h_phi_e[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(phi_e, w);
				h_Z[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(Z, w);
				h_Mx[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(M_x, w);
				h_Pt_pi[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(pT_pi, w);
				h_theta_pi[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(theta_pi, w);
				h_phi_pi[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(phi_pi, w);
				h_Vz_pi[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(Vz_pi, w);
				h_p_pi[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(p_pi, w);
				h_Phi_q[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(phi_q, w);
				h_Eta[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(eta, w);
				h_Xf[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill(xF, w);
				hQ2_Xb[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill( Q2, xB, w);
				hQ2_omega[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill( Q2, omega, w);
				hQ2_Z[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill( Q2, Z, w);
				hQ2_W[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill( xB, W, w);
				hBeta_p[chargeIdx][this_bin_Q2][this_bin_xB][v]->Fill( p_pi, beta, w);

				h_W[chargeIdx][this_bin_Q2][0][v]->Fill(W, w);
				h_Xb[chargeIdx][this_bin_Q2][0][v]->Fill(xB, w);
				h_Q2[chargeIdx][this_bin_Q2][0][v]->Fill(Q2, w);
				h_y[chargeIdx][this_bin_Q2][0][v]->Fill(y, w);
				h_omega[chargeIdx][this_bin_Q2][0][v]->Fill(omega, w);
				h_theta_e[chargeIdx][this_bin_Q2][0][v]->Fill(theta_e, w);
				h_Vz_e[chargeIdx][this_bin_Q2][0][v]->Fill(Vz_e, w);
				h_p_e[chargeIdx][this_bin_Q2][0][v]->Fill(p_e, w);
				h_phi_e[chargeIdx][this_bin_Q2][0][v]->Fill(phi_e, w);
				h_Z[chargeIdx][this_bin_Q2][0][v]->Fill(Z, w);
				h_Mx[chargeIdx][this_bin_Q2][0][v]->Fill(M_x, w);
				h_Pt_pi[chargeIdx][this_bin_Q2][0][v]->Fill(pT_pi, w);
				h_theta_pi[chargeIdx][this_bin_Q2][0][v]->Fill(theta_pi, w);
				h_phi_pi[chargeIdx][this_bin_Q2][0][v]->Fill(phi_pi, w);
				h_Vz_pi[chargeIdx][this_bin_Q2][0][v]->Fill(Vz_pi, w);
				h_p_pi[chargeIdx][this_bin_Q2][0][v]->Fill(p_pi, w);
				h_Phi_q[chargeIdx][this_bin_Q2][0][v]->Fill(phi_q, w);
				h_Eta[chargeIdx][this_bin_Q2][0][v]->Fill(eta, w);
				h_Xf[chargeIdx][this_bin_Q2][0][v]->Fill(xF, w);
				hQ2_Xb[chargeIdx][this_bin_Q2][0][v]->Fill( Q2, xB, w);
				hQ2_omega[chargeIdx][this_bin_Q2][0][v]->Fill( Q2, omega, w);
				hQ2_Z[chargeIdx][this_bin_Q2][0][v]->Fill( Q2, Z, w);
				hQ2_W[chargeIdx][this_bin_Q2][0][v]->Fill( xB, W, w);

				h_W[chargeIdx][0][this_bin_xB][v]->Fill(W, w);
				h_Xb[chargeIdx][0][this_bin_xB][v]->Fill(xB, w);
				h_Q2[chargeIdx][0][this_bin_xB][v]->Fill(Q2, w);
				h_y[chargeIdx][0][this_bin_xB][v]->Fill(y, w);
				h_omega[chargeIdx][0][this_bin_xB][v]->Fill(omega, w);
				h_theta_e[chargeIdx][0][this_bin_xB][v]->Fill(theta_e, w);
				h_Vz_e[chargeIdx][0][this_bin_xB][v]->Fill(Vz_e, w);
				h_p_e[chargeIdx][0][this_bin_xB][v]->Fill(p_e, w);
				h_phi_e[chargeIdx][0][this_bin_xB][v]->Fill(phi_e, w);
				h_Z[chargeIdx][0][this_bin_xB][v]->Fill(Z, w);
				h_Mx[chargeIdx][0][this_bin_xB][v]->Fill(M_x, w);
				h_Pt_pi[chargeIdx][0][this_bin_xB][v]->Fill(pT_pi, w);
				h_theta_pi[chargeIdx][0][this_bin_xB][v]->Fill(theta_pi, w);
				h_phi_pi[chargeIdx][0][this_bin_xB][v]->Fill(phi_pi, w);
				h_Vz_pi[chargeIdx][0][this_bin_xB][v]->Fill(Vz_pi, w);
				h_p_pi[chargeIdx][0][this_bin_xB][v]->Fill(p_pi, w);
				h_Phi_q[chargeIdx][0][this_bin_xB][v]->Fill(phi_q, w);
				h_Eta[chargeIdx][0][this_bin_xB][v]->Fill(eta, w);
				h_Xf[chargeIdx][0][this_bin_xB][v]->Fill(xF, w);
				hQ2_Xb[chargeIdx][0][this_bin_xB][v]->Fill( Q2, xB, w);
				hQ2_omega[chargeIdx][0][this_bin_xB][v]->Fill( Q2, omega, w);
				hQ2_Z[chargeIdx][0][this_bin_xB][v]->Fill( Q2, Z, w);
				hQ2_W[chargeIdx][0][this_bin_xB][v]->Fill( xB, W, w);
			}
		}
	}

	outFile->cd();

	for( int v = 0; v <= bins_var; v++ ){
		for( int j = 0; j <= bins_Q2; j++ ){
			for( int k = 0; k <= bins_xB; k++ ){
				for( int i = 0; i < 2; i++ ){
					h_W[i][j][k][v]->Write();
					h_Xb[i][j][k][v]->Write();
					h_Xf[i][j][k][v]->Write();
					h_Q2[i][j][k][v]->Write();
					h_Eta[i][j][k][v]->Write();
					h_y[i][j][k][v]->Write();
					h_omega[i][j][k][v]->Write();
					h_theta_e[i][j][k][v]->Write();
					h_Vz_e[i][j][k][v]->Write();
					h_p_e[i][j][k][v]->Write();
					h_phi_e[i][j][k][v]->Write();
					h_Z[i][j][k][v]->Write();
					h_p_pi[i][j][k][v]->Write();
					h_Vz_pi[i][j][k][v]->Write();
					h_theta_pi[i][j][k][v]->Write();
					h_phi_pi[i][j][k][v]->Write();
					h_Phi_q[i][j][k][v]->Write();
					h_Pt_pi[i][j][k][v]->Write();
					h_Mx[i][j][k][v]->Write();
					hQ2_Xb[i][j][k][v]->Write();
					hQ2_omega[i][j][k][v]->Write();
					hQ2_Z[i][j][k][v]->Write();
					hQ2_W[i][j][k][v]->Write();
					hBeta_p[i][j][k][v]->Write();
					hTheta_p[i][j][k][v]->Write();
				}
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
