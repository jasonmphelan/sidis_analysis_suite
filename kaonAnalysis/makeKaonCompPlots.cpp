#include <iostream>
#include <cmath>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "analyzer.h"
#include "reader.h"

using std::cerr;
using std::cout;
using namespace constants;
using namespace cutVals;

// kaon mass — note: Z and Mx in the skim tree are computed with Mpi (bug in test_skimmer),
// so those variables will be shifted. All other variables are unaffected.
static const double Mkaon = 0.493677; // GeV/c^2

int main( int argc, char** argv ){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 8 ){
		cerr << "Usage: ./makeKaonCompPlots [data_path] [mc_path] [output_name] "
		     << "[nFiles_data] [nFiles_mc] [EBeam] [target]\n";
		cerr << "  target: 0 = RGB/deuterium, 1 = RGA/proton\n";
		return -1;
	}

	TString data_path = argv[1];
	TString mc_path   = argv[2];
	TString out_name  = argv[3];
	int nFiles_data   = atoi(argv[4]);
	int nFiles_mc     = atoi(argv[5]);
	double EBeam      = atof(argv[6]);
	int target        = atoi(argv[7]);

	// -----------------------------------------------------------------------
	// Load files
	// -----------------------------------------------------------------------
	
	auto loadChain = [](TString path, int nFiles, int runType, double EBeam, int target) -> TChain* {
		TChain* chain = new TChain("ePi");
		reader r;
		r.setNumFiles(nFiles);
		r.setRunType(runType);
		r.setEnergy(EBeam);
		r.setTarget(target);
		r.getRunSkimsByName(chain, path);
		return chain;
	};

	TChain* chain_data = new TChain("ePi");//loadChain(data_path, nFiles_data, 0, EBeam, target);
	TChain* chain_mc   = new TChain("ePi");//loadChain(mc_path,   nFiles_mc,   1, EBeam, target);

	chain_data->Add("../trees/run_skim_6421.root");
	chain_mc->Add("../trees/final_skims/GEMC/detector_skim_10.2.root");

	// -----------------------------------------------------------------------
	// Analyzer (acceptance map, cuts)
	// -----------------------------------------------------------------------
	analyzer anal(0, -1);
	anal.setAnalyzerLevel(0);
	anal.setTarget(target);
	anal.loadAcceptanceMapContinuous(
		(TString)_DATA + "/acceptance_map/acceptanceMap_allE_final.root");

	// -----------------------------------------------------------------------
	// Histogram declarations
	// Index conventions:
	//   species: 0 = K+, 1 = K-
	//   source:  0 = data, 1 = MC
	// -----------------------------------------------------------------------
	const int NSPEC = 2;
	const int NSRC  = 2;
	TString spec_label[NSPEC] = {"kp", "km"};
	TString src_label[NSRC]   = {"data", "mc"};

	// Electron variables
	TH1F* h_Q2      [NSPEC][NSRC];
	TH1F* h_W       [NSPEC][NSRC];
	TH1F* h_xB      [NSPEC][NSRC];
	TH1F* h_y       [NSPEC][NSRC];
	TH1F* h_omega   [NSPEC][NSRC];
	TH1F* h_p_e     [NSPEC][NSRC];
	TH1F* h_theta_e [NSPEC][NSRC];
	TH1F* h_phi_e   [NSPEC][NSRC];
	TH1F* h_Vz_e    [NSPEC][NSRC];

	// Kaon kinematic variables
	TH1F* h_p_k      [NSPEC][NSRC];
	TH1F* h_theta_k  [NSPEC][NSRC];
	TH1F* h_phi_k    [NSPEC][NSRC];
	TH1F* h_pT_k     [NSPEC][NSRC];
	TH1F* h_Z_k      [NSPEC][NSRC]; // NOTE: computed with Mpi in skim — shifted
	TH1F* h_Mx_k     [NSPEC][NSRC]; // NOTE: same
	TH1F* h_eta_k    [NSPEC][NSRC];
	TH1F* h_Vz_k     [NSPEC][NSRC];
	TH1F* h_dVz_k    [NSPEC][NSRC]; // Vz_k - Vz_e

	// Kaon detector variables
	TH1F* h_Chi2_k   [NSPEC][NSRC];
	TH1F* h_Beta_k   [NSPEC][NSRC]; // TOF beta
	TH1F* h_m2_k     [NSPEC][NSRC]; // TOF mass^2
	TH1F* h_Beta_rich[NSPEC][NSRC];
	TH1F* h_m2_rich  [NSPEC][NSRC]; // RICH mass^2
	TH1F* h_edge0    [NSPEC][NSRC];
	TH1F* h_edge1    [NSPEC][NSRC];
	TH1F* h_edge2    [NSPEC][NSRC];
	TH1F* h_DC_chi2  [NSPEC][NSRC];
	TH1F* h_DC_sec   [NSPEC][NSRC];

	// 2D
	TH2F* h_Q2_W     [NSPEC][NSRC];
	TH2F* h_m2_p     [NSPEC][NSRC]; // TOF m^2 vs p
	TH2F* h_m2rich_p [NSPEC][NSRC]; // RICH m^2 vs p
	TH2F* h_Beta_p   [NSPEC][NSRC];
	TH2F* h_theta_p  [NSPEC][NSRC]; // theta vs p

	for( int s = 0; s < NSPEC; s++ ){
		for( int src = 0; src < NSRC; src++ ){
			TString tag = spec_label[s] + "_" + src_label[src];

			// Electron
			h_Q2[s][src]      = new TH1F("hQ2_"     +tag, ";Q^{2} [GeV^{2}];Counts [a.u.]",  60, 2, 8);
			h_W[s][src]       = new TH1F("hW_"      +tag, ";W [GeV];Counts [a.u.]",           60, 2, 4.5);
			h_xB[s][src]      = new TH1F("hxB_"     +tag, ";x_{B};Counts [a.u.]",             50, 0.1, 0.66);
			h_y[s][src]       = new TH1F("hY_"      +tag, ";y;Counts [a.u.]",                 50, 0, 0.75);
			h_omega[s][src]   = new TH1F("hOmega_"  +tag, ";#nu [GeV];Counts [a.u.]",         50, 0, 10);
			h_p_e[s][src]     = new TH1F("hP_e_"    +tag, ";p_{e} [GeV/c];Counts [a.u.]",     50, 3, 11);
			h_theta_e[s][src] = new TH1F("hTheta_e_"+tag, ";#theta_{e} [deg];Counts [a.u.]",  60, 5, 35);
			h_phi_e[s][src]   = new TH1F("hPhi_e_"  +tag, ";#phi_{e} [deg];Counts [a.u.]",   120,-180,180);
			h_Vz_e[s][src]    = new TH1F("hVz_e_"   +tag, ";V_{z}^{e} [cm];Counts [a.u.]",   60,-15,5);

			// Kaon kinematics
			h_p_k[s][src]     = new TH1F("hP_k_"    +tag, ";p_{K} [GeV/c];Counts [a.u.]",    60, 1.25, 5);
			h_theta_k[s][src] = new TH1F("hTheta_k_"+tag, ";#theta_{K} [deg];Counts [a.u.]",  60, 5, 35);
			h_phi_k[s][src]   = new TH1F("hPhi_k_"  +tag, ";#phi_{K} [deg];Counts [a.u.]",   120,-180,180);
			h_pT_k[s][src]    = new TH1F("hPt_k_"   +tag, ";p_{T,K} [GeV/c];Counts [a.u.]",  50, 0, 1.3);
			h_Z_k[s][src]     = new TH1F("hZ_k_"    +tag, ";z_{K} (Mpi mass);Counts [a.u.]",  50, 0.3, 1.0);
			h_Mx_k[s][src]    = new TH1F("hMx_k_"   +tag, ";M_{X} [GeV] (Mpi mass);Counts [a.u.]", 60, 1.5, 4);
			h_eta_k[s][src]   = new TH1F("hEta_k_"  +tag, ";#eta_{K};Counts [a.u.]",          50, 0, 4.5);
			h_Vz_k[s][src]    = new TH1F("hVz_k_"   +tag, ";V_{z}^{K} [cm];Counts [a.u.]",   60,-15,5);
			h_dVz_k[s][src]   = new TH1F("hdVz_k_"  +tag, ";#DeltaV_{z} [cm];Counts [a.u.]", 60,-15,10);

			// Kaon detector
			h_Chi2_k[s][src]    = new TH1F("hChi2_k_"   +tag, ";#chi^{2}_{PID};Counts [a.u.]",   80,-5,5);
			h_Beta_k[s][src]    = new TH1F("hBeta_k_"   +tag, ";#beta_{TOF};Counts [a.u.]",       80, 0.7,1.05);
			h_m2_k[s][src]      = new TH1F("hm2_k_"     +tag, ";m^{2}_{TOF} [GeV^{2}];Counts [a.u.]", 100,-0.2, 0.8);
			h_Beta_rich[s][src] = new TH1F("hBeta_rich_"+tag, ";#beta_{RICH};Counts [a.u.]",      80, 0.7,1.05);
			h_m2_rich[s][src]   = new TH1F("hm2_rich_"  +tag, ";m^{2}_{RICH} [GeV^{2}];Counts [a.u.]",100,-0.2,0.8);
			h_edge0[s][src]     = new TH1F("hEdge0_k_"  +tag, ";DC edge R1 [cm];Counts [a.u.]",   50, 0, 30);
			h_edge1[s][src]     = new TH1F("hEdge1_k_"  +tag, ";DC edge R2 [cm];Counts [a.u.]",   50, 0, 30);
			h_edge2[s][src]     = new TH1F("hEdge2_k_"  +tag, ";DC edge R3 [cm];Counts [a.u.]",   50, 0, 30);
			h_DC_chi2[s][src]   = new TH1F("hDC_chi2_k_"+tag, ";DC #chi^{2}/NDF;Counts [a.u.]",   60, 0, 30);
			h_DC_sec[s][src]    = new TH1F("hDC_sec_k_" +tag, ";DC sector;Counts [a.u.]",          7, 0, 7);

			// 2D
			h_Q2_W[s][src]    = new TH2F("hQ2_W_"    +tag, ";Q^{2} [GeV^{2}];W [GeV]",         50,2,8, 50,2,4.5);
			h_m2_p[s][src]    = new TH2F("hm2_p_"    +tag, ";p_{K} [GeV/c];m^{2}_{TOF}",        50,1.25,5, 100,-0.2,0.8);
			h_m2rich_p[s][src]= new TH2F("hm2rich_p_"+tag, ";p_{K} [GeV/c];m^{2}_{RICH}",       50,1.25,5, 100,-0.2,0.8);
			h_Beta_p[s][src]  = new TH2F("hBeta_p_"  +tag, ";p_{K} [GeV/c];#beta_{TOF}",        50,1.25,5,  80,0.7,1.05);
			h_theta_p[s][src] = new TH2F("hTheta_p_" +tag, ";p_{K} [GeV/c];#theta_{K} [deg]",   50,1.25,5,  60,5,35);
		}
	}

	// -----------------------------------------------------------------------
	// Event loop — run over data (src=0) then MC (src=1)
	// -----------------------------------------------------------------------
	TChain* chains[NSRC] = { chain_data, chain_mc };

	for( int src = 0; src < NSRC; src++ ){

		TTreeReader treader( chains[src] );
		TTreeReaderValue<electron> e( treader, "e" );
		TString branch = (src == 0) ? "ka" : "pi";
		TTreeReaderArray<pion>     k( treader, branch );

		long long event_count = 0;
		long long event_total = treader.GetEntries();
		cout << "Source " << src_label[src] << ": " << event_total << " events\n";

		while( treader.Next() ){
			if( event_count % 500000 == 0 ){
				cout << src_label[src] << ": " << event_count << " / " << event_total << "\n";
			}
			event_count++;

			// Electron acceptance map
			double p_e     = e->get3Momentum().Mag();
			double theta_e = e->get3Momentum().Theta() * rad_to_deg;
			double phi_e   = e->get3Momentum().Phi()   * rad_to_deg;
			if( anal.applyAcceptanceMap(p_e, phi_e, theta_e, 0) < 0 ) continue;

			// Basic DIS cuts
			double Q2 = e->getQ2();
			double W  = sqrt( e->getW2() );
			double y  = e->getY();
			if( W  < W_min  ) continue;
			if( Q2 < Q2_min ) continue;
			if( y  > y_max  ) continue;

			int nKaons = (int)( k.end() - k.begin() );
			for( int i = 0; i < nKaons; i++ ){

				if( src == 1 && abs(k[i].getPID_eb()) != 321 ) continue;
			//if( src == 1 && !anal.applyPionDetectorChi2(k[i]) ) continue;

				int s = (int)( k[i].getCharge() < 0 ); // 0 = K+, 1 = K-

				double p_k     = k[i].get3Momentum().Mag();
				double theta_k = k[i].get3Momentum().Theta() * rad_to_deg;
				double phi_k   = k[i].get3Momentum().Phi()   * rad_to_deg;

				// Momentum and angular acceptance
				if( p_k < P_pi_min || p_k > P_pi_max ) continue;
				if( theta_k < theta_min || theta_k > theta_max ) continue;
				if( anal.applyAcceptanceMap(p_k, phi_k, theta_k, s+1) < 0 ) continue;

				// Mx cut (note: stored with Mpi mass assumption)
				//if( k[i].getMx() < Mx_min ) continue;

				// Compute derived quantities
				double beta_tof  = k[i].getBeta();
				double beta_rich = k[i].getBeta_rich();

				double m2_tof  = (beta_tof  > 0.001) ?
					p_k*p_k * (1 - beta_tof *beta_tof ) / (beta_tof *beta_tof ) : -999;
				double m2_rich = (beta_rich > 0.001) ?
					p_k*p_k * (1 - beta_rich*beta_rich) / (beta_rich*beta_rich) : -999;

				double dVz = ( k[i].getVt() - e->getVt() ).Z();

				// ---- fill electron variables ----
				h_Q2[s][src]     ->Fill( Q2 );
				h_W[s][src]      ->Fill( W );
				h_xB[s][src]     ->Fill( e->getXb() );
				h_y[s][src]      ->Fill( y );
				h_omega[s][src]  ->Fill( e->getOmega() );
				h_p_e[s][src]    ->Fill( p_e );
				h_theta_e[s][src]->Fill( theta_e );
				h_phi_e[s][src]  ->Fill( phi_e );
				h_Vz_e[s][src]   ->Fill( e->getVt().Z() );

				h_Q2_W[s][src]   ->Fill( Q2, W );
				h_theta_p[s][src]->Fill( p_k, theta_k );

				// ---- fill kaon kinematic variables ----
				h_p_k[s][src]    ->Fill( p_k );
				h_theta_k[s][src]->Fill( theta_k );
				h_phi_k[s][src]  ->Fill( phi_k );
				h_pT_k[s][src]   ->Fill( k[i].getPi_q().Pt() );
				h_Z_k[s][src]    ->Fill( k[i].getZ() );
				h_Mx_k[s][src]   ->Fill( k[i].getMx() );
				h_eta_k[s][src]  ->Fill( k[i].getEta() );
				h_Vz_k[s][src]   ->Fill( k[i].getVt().Z() );
				h_dVz_k[s][src]  ->Fill( dVz );

				// ---- fill kaon detector variables ----
				h_Chi2_k[s][src]   ->Fill( k[i].getChi2() );
				h_Beta_k[s][src]   ->Fill( beta_tof );
				h_Beta_p[s][src]   ->Fill( p_k, beta_tof );
				if( m2_tof  > -998 ) h_m2_k[s][src]   ->Fill( m2_tof );
				if( m2_tof  > -998 ) h_m2_p[s][src]   ->Fill( p_k, m2_tof );
				if( beta_rich > 0.001 ){
					h_Beta_rich[s][src]->Fill( beta_rich );
					if( m2_rich > -998 ) h_m2_rich[s][src] ->Fill( m2_rich );
					if( m2_rich > -998 ) h_m2rich_p[s][src]->Fill( p_k, m2_rich );
				}
				h_edge0[s][src]  ->Fill( k[i].getEdge(0) );
				h_edge1[s][src]  ->Fill( k[i].getEdge(1) );
				h_edge2[s][src]  ->Fill( k[i].getEdge(2) );
				if( k[i].getDC_NDF() > 0 ){
					h_DC_chi2[s][src]->Fill( k[i].getDC_chi2() / k[i].getDC_NDF() );
				}
				h_DC_sec[s][src] ->Fill( k[i].getDC_sector() );
			}
		}
	} // end source loop

	// -----------------------------------------------------------------------
	// Style
	// -----------------------------------------------------------------------
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);

	// -----------------------------------------------------------------------
	// Normalize MC to data area, write histograms, draw canvases
	// -----------------------------------------------------------------------
	TFile* outFile = new TFile( out_name + ".root", "RECREATE" );
	outFile->cd();

	// Style helpers
	auto styleData = []( TH1F* h ){
		h->SetMarkerStyle(20);
		h->SetMarkerSize(0.7);
		h->SetMarkerColor(kBlack);
		h->SetLineColor(kBlack);
	};
	auto styleMC = []( TH1F* h ){
		h->SetLineColor(kRed+1);
		h->SetFillColorAlpha(kRed+1, 0.25);
		h->SetLineWidth(2);
	};

	// Each entry: { hist array[NSPEC][NSRC], canvas name, group label }
	struct VarEntry {
		TH1F* (*arr)[NSRC];
		TString name;
	};

	// Groups of related variables — each group goes on one canvas (2 rows x N cols)
	// Row 0 = K+, Row 1 = K-
	struct Group {
		TString title;
		std::vector<VarEntry> vars;
	};

	std::vector<Group> groups = {
		{ "Electron kinematics",
		  { {h_Q2, "Q^{2}"}, {h_W, "W"}, {h_xB, "x_{B}"},
		    {h_y, "y"}, {h_omega, "#nu"} } },
		{ "Electron tracking",
		  { {h_p_e, "p_{e}"}, {h_theta_e, "#theta_{e}"},
		    {h_phi_e, "#phi_{e}"}, {h_Vz_e, "V_{z}^{e}"} } },
		{ "Kaon kinematics",
		  { {h_p_k, "p_{K}"}, {h_theta_k, "#theta_{K}"},
		    {h_phi_k, "#phi_{K}"}, {h_pT_k, "p_{T,K}"},
		    {h_Z_k, "z_{K}"}, {h_Mx_k, "M_{X}"}, {h_eta_k, "#eta_{K}"} } },
		{ "Kaon vertex",
		  { {h_Vz_k, "V_{z}^{K}"}, {h_dVz_k, "#DeltaV_{z}"} } },
		{ "Kaon PID",
		  { {h_Chi2_k, "#chi^{2}_{PID}"}, {h_Beta_k, "#beta_{TOF}"},
		    {h_m2_k, "m^{2}_{TOF}"}, {h_Beta_rich, "#beta_{RICH}"},
		    {h_m2_rich, "m^{2}_{RICH}"} } },
		{ "Kaon DC fiducials",
		  { {h_edge0, "Edge R1"}, {h_edge1, "Edge R2"},
		    {h_edge2, "Edge R3"}, {h_DC_chi2, "DC #chi^{2}/NDF"},
		    {h_DC_sec, "DC sector"} } },
	};

	// Normalize all MC histograms to data before drawing
	for( int s = 0; s < NSPEC; s++ ){
		TH1F** all_pairs[] = {
			h_Q2[s], h_W[s], h_xB[s], h_y[s], h_omega[s],
			h_p_e[s], h_theta_e[s], h_phi_e[s], h_Vz_e[s],
			h_p_k[s], h_theta_k[s], h_phi_k[s], h_pT_k[s],
			h_Z_k[s], h_Mx_k[s], h_eta_k[s], h_Vz_k[s], h_dVz_k[s],
			h_Chi2_k[s], h_Beta_k[s], h_m2_k[s],
			h_Beta_rich[s], h_m2_rich[s],
			h_edge0[s], h_edge1[s], h_edge2[s],
			h_DC_chi2[s], h_DC_sec[s]
		};
		for( auto& pair : all_pairs ){
			styleData( pair[0] );
			styleMC(   pair[1] );
			if( pair[0]->Integral() > 0 && pair[1]->Integral() > 0 )
				pair[1]->Scale( pair[0]->Integral() / pair[1]->Integral() );
			pair[0]->Write();
			pair[1]->Write();
		}
		for( int src = 0; src < NSRC; src++ ){
			h_Q2_W[s][src]   ->Write();
			h_m2_p[s][src]   ->Write();
			h_m2rich_p[s][src]->Write();
			h_Beta_p[s][src] ->Write();
			h_theta_p[s][src]->Write();
		}
	}

	// Draw canvases: one per group, 2 rows (K+/K-) x N cols (variables)
	TString charge_label[NSPEC] = {"K^{+}", "K^{-}"};

	for( auto& grp : groups ){
		int nVars = (int)grp.vars.size();
		int nCols = nVars;
		int nRows = NSPEC; // 2 rows: K+ and K-

		int cW = std::max(800, nCols * 300);
		int cH = nRows * 300;

		TCanvas* c = new TCanvas( "c_" + grp.title, grp.title, cW, cH );
		c->Divide( nCols, nRows );

		for( int s = 0; s < NSPEC; s++ ){
			for( int v = 0; v < nVars; v++ ){
				int pad = s * nCols + v + 1;
				c->cd(pad);
				gPad->SetLeftMargin(0.15);
				gPad->SetBottomMargin(0.15);

				TH1F* hd = grp.vars[v].arr[s][0]; // data
				TH1F* hm = grp.vars[v].arr[s][1]; // MC

				// Set y-axis range to accommodate both
				double ymax = std::max( hd->GetMaximum(), hm->GetMaximum() ) * 1.3;
				hd->GetYaxis()->SetRangeUser(0, ymax);
				hd->GetYaxis()->SetTitle("Counts [a.u.]");
				hd->GetYaxis()->SetTitleSize(0.055);
				hd->GetYaxis()->SetLabelSize(0.050);
				hd->GetXaxis()->SetTitle( grp.vars[v].name );
				hd->GetXaxis()->SetTitleSize(0.055);
				hd->GetXaxis()->SetLabelSize(0.050);

				hd->Draw("E");
				hm->Draw("HIST SAME");
				hd->Draw("E SAME"); // data on top

				// Column header (variable name) in top row only
				if( s == 0 ){
					TLatex* lat = new TLatex();
					lat->SetNDC();
					lat->SetTextSize(0.06);
					lat->SetTextAlign(22);
					lat->DrawLatex(0.5, 0.96, grp.vars[v].name);
				}

				// Row label (charge) in leftmost column only
				if( v == 0 ){
					TLatex* lat = new TLatex();
					lat->SetNDC();
					lat->SetTextSize(0.06);
					lat->SetTextAlign(22);
					lat->SetTextAngle(90);
					lat->DrawLatex(0.04, 0.5, charge_label[s]);
				}

				// Legend in top-left variable of each row
				if( v == 0 ){
					TLegend* leg = new TLegend(0.55, 0.72, 0.92, 0.90);
					leg->SetBorderSize(0);
					leg->SetFillStyle(0);
					leg->SetTextSize(0.052);
					leg->AddEntry(hd, "Data", "lep");
					leg->AddEntry(hm, "MC",   "lf");
					leg->Draw();
				}
			}
		}

		c->Write();
		c->SaveAs( out_name + "_" + grp.title + ".pdf" );
		delete c;
	}

	// Theta-p canvas: 2 rows (K+/K-) x 2 cols (data/MC), COLZ
	{
		TCanvas* c = new TCanvas("c_theta_p", "#theta vs p", 1000, 700);
		c->Divide(2, NSPEC);
		TString src_title[NSRC] = {"Data", "MC"};
		for( int s = 0; s < NSPEC; s++ ){
			for( int src = 0; src < NSRC; src++ ){
				c->cd( s*2 + src + 1 );
				gPad->SetLeftMargin(0.14);
				gPad->SetRightMargin(0.14);
				gPad->SetBottomMargin(0.14);
				h_theta_p[s][src]->GetXaxis()->SetTitleSize(0.055);
				h_theta_p[s][src]->GetYaxis()->SetTitleSize(0.055);
				h_theta_p[s][src]->GetXaxis()->SetLabelSize(0.050);
				h_theta_p[s][src]->GetYaxis()->SetLabelSize(0.050);
				h_theta_p[s][src]->Draw("COLZ");
				TLatex lat;
				lat.SetNDC();
				lat.SetTextSize(0.055);
				lat.SetTextAlign(22);
				lat.DrawLatex(0.5, 0.96,
					charge_label[s] + " " + src_title[src]);
			}
		}
		c->Write();
		c->SaveAs( out_name + "_theta_p.pdf" );
		delete c;
	}

	outFile->Close();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Done. Elapsed time: " << elapsed.count() << " s\n";
	return 0;
}
