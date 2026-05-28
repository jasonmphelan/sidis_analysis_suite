#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "analyzer.h"
#define HIST_PATH _HIST

using std::cerr;
using std::isfinite;
using std::cout;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv ){

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		cerr << "[nPhotons] [prob threshold (optional)] \n";
		return -1;
	}
	cerr << "Output file: " << argv[1] << "\n";

	TString out_name  = argv[1];
	int    nPhotons   = atoi(argv[2]);
	double prob       = (argc > 3) ? atof(argv[3]) : 0.0;

	// ── Parallelization setup ────────────────────────────────────────────
	ROOT::EnableImplicitMT();

	// ── Input files ──────────────────────────────────────────────────────
	std::vector<std::string> files = {
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6421.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6422.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6426.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6429.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6430.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6431.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6432.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6433.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6437.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6442.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6444.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6445.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6449.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6450.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6451.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6452.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6453.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6454.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6455.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6456.root",
		"/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_6457.root",
	};

	ROOT::RDataFrame df_full("ePi", files);
	// Mirror the 100M event cap from the original
	auto df = df_full.Range(100000000ULL);

	int nSlots = df.GetNSlots();
	cout << "Running with " << nSlots << " thread slots\n";

	// ── Per-slot histograms ──────────────────────────────────────────────
	int    p_bins = 20;
	double p_lo   = 1.25, p_hi = 5.0;
	TString species[4] = {"pip", "pim", "kp", "km"};

	// slot → [4 species]
	std::vector<std::array<TH1F*, 4>> hP_eb_s    (nSlots);
	std::vector<std::array<TH1F*, 4>> hP_true_s  (nSlots);
	std::vector<std::array<TH1F*, 4>> hP_opp_s   (nSlots);
	std::vector<std::array<TH2F*, 4>> hTheta_p_eb_s   (nSlots);
	std::vector<std::array<TH2F*, 4>> hTheta_p_true_s (nSlots);

	for( int s = 0; s < nSlots; s++ ){
		for( int i = 0; i < 4; i++ ){
			TString tag = Form("_s%d_", s);
			hP_eb_s[s][i] = new TH1F(
				"hP_eb_"   + species[i] + tag,
				"p (" + species[i] + ");p (GeV/c);Counts",
				p_bins, p_lo, p_hi );
			hP_true_s[s][i] = new TH1F(
				"hP_true_" + species[i] + tag,
				"p (" + species[i] + ");p (GeV/c);Counts",
				p_bins, p_lo, p_hi );
			hP_opp_s[s][i] = new TH1F(
				"hP_opp_"  + species[i] + tag,
				"p (" + species[i] + ");p (GeV/c);Counts",
				p_bins, p_lo, p_hi );
			hTheta_p_eb_s[s][i] = new TH2F(
				"hTheta_p_eb_"   + species[i] + tag,
				"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
				100, p_lo, p_hi, 75, 5.0, 35.0 );
			hTheta_p_true_s[s][i] = new TH2F(
				"hTheta_p_true_" + species[i] + tag,
				"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
				100, p_lo, p_hi, 75, 5.0, 35.0 );
		}
	}

	// ── Per-slot analyzers ───────────────────────────────────────────────
	// Each slot needs its own analyzer because TF1::Eval is not thread-safe
	// when called on shared objects (internal caching).
	std::vector<analyzer*> anal_slots(nSlots);
	for( int s = 0; s < nSlots; s++ ){
		anal_slots[s] = new analyzer(0, -1);
		anal_slots[s]->setAnalyzerLevel(0);
		anal_slots[s]->loadAcceptanceMapContinuous(
			(TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root" );
	}

	// ── Parallel event loop via ForeachSlot ──────────────────────────────
	cout << "Beginning Event Loop\n";

	df.ForeachSlot(
		[&]( unsigned int slot,
		     const electron&            e,
		     const ROOT::RVec<pion>&    pi,
		     const ROOT::RVec<pion>&    k,
		     const ROOT::RVec<int>&     pidRICH_pi,
		     const ROOT::RVec<int>&     pidRICH_k,
		     const ROOT::RVec<double>&  probRICH_pi,
		     const ROOT::RVec<double>&  probRICH_k,
		     const ROOT::RVec<double>&  nRICH_pi,
		     const ROOT::RVec<double>&  nRICH_k )
		{
			auto& anal = *anal_slots[slot];

			double p_e     = e.get3Momentum().Mag();
			double theta_e = e.get3Momentum().Theta() * rad_to_deg;
			double phi_e   = e.get3Momentum().Phi()   * rad_to_deg;

			if( anal.applyAcceptanceMap(p_e, phi_e, theta_e, 0) < 0 ) return;

			// ── Pion loop ────────────────────────────────────────────
			for( int i = 0; i < (int)pi.size(); i++ ){
				if( pi[i].getBeta_rich() < .0001 )       continue;
				if( nRICH_pi[i] <= nPhotons - 1 )        continue;
				if( probRICH_pi[i] < prob )               continue;

				int    chargeIdx = (int)(pi[i].getCharge() < 1);
				double p_pi      = pi[i].get3Momentum().Mag();
				double theta_pi  = pi[i].get3Momentum().Theta() * rad_to_deg;

				if( p_pi < p_lo || p_pi > p_hi ) continue;

				if( abs(pidRICH_pi[i]) == 211 ){
					hP_true_s[slot][chargeIdx]->Fill(p_pi);
					hTheta_p_true_s[slot][chargeIdx]->Fill(p_pi, theta_pi);
					hP_eb_s[slot][chargeIdx]->Fill(p_pi);
					hTheta_p_eb_s[slot][chargeIdx]->Fill(p_pi, theta_pi);
				}
				if( abs(pidRICH_pi[i]) == 321 ){
					hP_true_s[slot][2 + chargeIdx]->Fill(p_pi);
					hP_opp_s[slot][2 + chargeIdx]->Fill(p_pi);
					hTheta_p_true_s[slot][2 + chargeIdx]->Fill(p_pi, theta_pi);
				}
			}

			// ── Kaon loop ────────────────────────────────────────────
			for( int i = 0; i < (int)k.size(); i++ ){
				if( k[i].getBeta_rich() < .0001 )        continue;
				if( nRICH_k[i] <= nPhotons - 1 )         continue;
				if( probRICH_k[i] < prob )                continue;

				int    chargeIdx = (int)(k[i].getCharge() < 1);
				double p_k       = k[i].get3Momentum().Mag();
				double theta_k   = k[i].get3Momentum().Theta() * rad_to_deg;

				if( p_k < p_lo || p_k > p_hi ) continue;

				if( abs(pidRICH_k[i]) == 321 ){
					hP_true_s[slot][chargeIdx + 2]->Fill(p_k);
					hTheta_p_true_s[slot][chargeIdx + 2]->Fill(p_k, theta_k);
					hP_eb_s[slot][chargeIdx + 2]->Fill(p_k);
					hTheta_p_eb_s[slot][chargeIdx + 2]->Fill(p_k, theta_k);
				}
				if( abs(pidRICH_k[i]) == 211 ){
					hP_opp_s[slot][chargeIdx]->Fill(p_k);
					hP_true_s[slot][chargeIdx]->Fill(p_k);
					hTheta_p_true_s[slot][chargeIdx]->Fill(p_k, theta_k);
				}
			}
		},
		{"e", "pi", "ka", "pidRICH", "pidRICH_ka",
		 "probRICH", "probRICH_ka", "nRICH", "nRICH_ka"}
	);

	// ── Merge per-slot histograms into slot-0 ────────────────────────────
	for( int s = 1; s < nSlots; s++ ){
		for( int i = 0; i < 4; i++ ){
			hP_eb_s[0][i]        ->Add( hP_eb_s[s][i]         );
			hP_true_s[0][i]      ->Add( hP_true_s[s][i]       );
			hP_opp_s[0][i]       ->Add( hP_opp_s[s][i]        );
			hTheta_p_eb_s[0][i]  ->Add( hTheta_p_eb_s[s][i]   );
			hTheta_p_true_s[0][i]->Add( hTheta_p_true_s[s][i] );
		}
	}

	// ── Write output ─────────────────────────────────────────────────────
	TFile* outFile = new TFile(out_name, "RECREATE");
	outFile->cd();

	for( int i = 0; i < 4; i++ ){
		hP_true_s[0][i]      ->SetName( "hP_true_"      + species[i] ); hP_true_s[0][i]      ->Write();
		hP_eb_s[0][i]        ->SetName( "hP_eb_"        + species[i] ); hP_eb_s[0][i]        ->Write();
		hTheta_p_eb_s[0][i]  ->SetName( "hTheta_p_eb_"  + species[i] ); hTheta_p_eb_s[0][i]  ->Write();
		hTheta_p_true_s[0][i]->SetName( "hTheta_p_true_"+ species[i] ); hTheta_p_true_s[0][i]->Write();
		hP_opp_s[0][i]       ->SetName( "hP_opp_"       + species[i] ); hP_opp_s[0][i]       ->Write();
	}

	TString ratio_species[4] = {"pip", "pim", "kp", "km"};
	for( int i = 0; i < 4; i++ ){
		TH1F* hRatio_p = (TH1F*)hP_eb_s[0][i]->Clone("hRatio_p_" + ratio_species[i]);
		hRatio_p->Divide( hP_true_s[0][i] );
		hRatio_p->SetTitle("p^{+}/p^{-} ratio (" + ratio_species[i] + ");p (GeV/c);N^{+}/N^{-}");
		hRatio_p->Write();

		TH2F* hRatio_theta_p = (TH2F*)hTheta_p_eb_s[0][i]->Clone("hRatio_theta_p_" + ratio_species[i]);
		hRatio_theta_p->Divide( hTheta_p_true_s[0][i] );
		hRatio_theta_p->SetTitle("#theta vs p ratio (" + ratio_species[i] + ");p (GeV/c);#theta (deg)");
		hRatio_theta_p->Write();
	}

	cout << "FINISHED WRITING\n";
	outFile->Close();

	for( int s = 0; s < nSlots; s++ ) delete anal_slots[s];

	return 0;
}
