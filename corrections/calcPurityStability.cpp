#include <fstream>
#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH3.h>
#include <TChain.h>
#include "electron.h"
#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
#include "constants.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

// Purity  = N_match(rec bin == gen bin) / N_rec(bin)
// Stability = N_match(rec bin == gen bin) / N_gen(bin)
//
// Usage: ./calcPurityStability [Rec File] [Gen File] [Output File]
//                              [Matching] [Target] [theta_cut]
//                              (opt: matching file)
//
// Rec File: output of final_skimmer (has both pi and pi_gen branches)
// Gen File: output of gen_skimmer   (has only pi_gen branch)

using std::cerr;
using std::cout;
using namespace cutVals;
using namespace constants;

const int nXbBins = bins_xB;
const int nQ2Bins = bins_Q2;
const int nZBins  = bins_Z;

double getUncertainty( double num, double den );

int main( int argc, char** argv ){

	if( argc < 7 ){
		cerr << "Usage: ./calcPurityStability [Rec File] [Gen File] [Output File]\n";
		cerr << "                             [Matching] [Target] [theta_cut]\n";
		cerr << "                             (opt: matching file)\n";
		return -1;
	}

	TString inName_rec   = argv[1];
	TString inName_gen   = argv[2];
	TString outName      = argv[3];
	int     matchType    = atoi( argv[4] );
	int     rga          = atoi( argv[5] );
	double  theta_cut    = atof( argv[6] );
	TString matchingFile = "matchCut2D_map.root";
	if( argc > 7 ) matchingFile = argv[7];

	cerr << "Rec: " << inName_rec << "\nGen: " << inName_gen << "\n";

	TFile* outFile    = new TFile( outName, "RECREATE" );
	TFile* inFile_rec = new TFile( inName_rec );
	TTree* recChain   = (TTree*) inFile_rec->Get("ePi");
	TFile* inFile_gen = new TFile( inName_gen );
	TTree* genChain   = (TTree*) inFile_gen->Get("ePi");

	// ── Output histograms ────────────────────────────────────────────────────
	// N_rec:   reconstructed events per (xB, Q2, z) bin
	// N_gen:   generated events per (xB, Q2, z) bin
	// N_match: events where rec bin == gen bin (same xB, Q2, z)
	TH3F* hRec_p   = new TH3F("hRec_p",   "hRec_p",   nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hRec_m   = new TH3F("hRec_m",   "hRec_m",   nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hGen_p   = new TH3F("hGen_p",   "hGen_p",   nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hGen_m   = new TH3F("hGen_m",   "hGen_m",   nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hMatch_p = new TH3F("hMatch_p", "hMatch_p", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hMatch_m = new TH3F("hMatch_m", "hMatch_m", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);

	TH3F* hPurity_p    = new TH3F("hPurity_p",    "Purity (#pi^{+});x_{B};Q^{2};z",    nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hPurity_m    = new TH3F("hPurity_m",    "Purity (#pi^{-});x_{B};Q^{2};z",    nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hStability_p = new TH3F("hStability_p", "Stability (#pi^{+});x_{B};Q^{2};z", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);
	TH3F* hStability_m = new TH3F("hStability_m", "Stability (#pi^{-});x_{B};Q^{2};z", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1.);

	// ── Analyzer ─────────────────────────────────────────────────────────────
	analyzer anal(0, -1);
	anal.setAnalyzerLevel(1);
	anal.loadMatchingFunctions(matchingFile);
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous(
		(TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");

	// ── Reconstructed tree loop ───────────────────────────────────────────────
	// Fills hRec (all selected rec pions) and hMatch (rec pions whose gen
	// counterpart lands in the same xB, Q2, z bin).
	TTreeReader reader_rec(recChain);
	TTreeReaderArray<bool>      isGoodPion(        reader_rec, "isGoodPion");
	TTreeReaderArray<bool>      isGoodPion_no_acc( reader_rec, "isGoodPion_no_acc");
	TTreeReaderArray<bool>      isGoodGenPion(     reader_rec, "isGoodGenPion");
	TTreeReaderArray<bool>      isGoodPion3d(      reader_rec, "isGoodPion_3d");
	TTreeReaderValue<electron>  e(     reader_rec, "e");
	TTreeReaderValue<genElectron> e_MC(reader_rec, "e_gen");
	TTreeReaderArray<pion>      pi_vec(  reader_rec, "pi");
	TTreeReaderArray<genPion>   pi_match(reader_rec, "pi_gen");

	int event_total = recChain->GetEntries();
	cout << "Processing rec tree (" << event_total << " events)...\n";

	while( reader_rec.Next() ){
		int event_count = reader_rec.GetCurrentEntry();
		if( event_count % 100000 == 0 )
			cout << "  " << event_count << " / " << event_total << "\n";

		double Q2_rec = e->getQ2();
		double xB_rec = e->getXb();
		double Q2_gen = e_MC->getQ2();
		double xB_gen = e_MC->getXb();

		int bQ2_rec = (int)( (Q2_rec - Q2_min) / (Q2_max - Q2_min) * nQ2Bins );
		int bxB_rec = (int)( (xB_rec - xB_min) / (xB_max - xB_min) * nXbBins );
		int bQ2_gen = (int)( (Q2_gen - Q2_min) / (Q2_max - Q2_min) * nQ2Bins );
		int bxB_gen = (int)( (xB_gen - xB_min) / (xB_max - xB_min) * nXbBins );

		int pi_count = -1;
		for( auto pi : pi_vec ){
			pi_count++;

			if( abs( pi.getPID() ) != 211 ) continue;

			int    chargeIdx = (int)( pi.getCharge() < 1 );
			double p_mom     = pi.get3Momentum().Mag();
			double p_theta   = pi.get3Momentum().Theta() * rad_to_deg;

			if( p_theta < theta_cut ) continue;

			bool matching = true;
			if( matchType == 1 || matchType == 2 )
				matching = anal.applyAcceptanceMatching(pi, 2);
			else if( matchType == 3 )
				matching = anal.applyAcceptanceMap( p_mom, rad_to_deg*pi.get3Momentum().Phi(), p_theta, 1 ) >= 0 &&
				           anal.applyAcceptanceMap( p_mom, rad_to_deg*pi.get3Momentum().Phi(), p_theta, 2 ) >= 0;

			bool selected = ( matchType==0 && isGoodPion_no_acc[pi_count] ) ||
			                ( matchType==1 && matching && isGoodPion_no_acc[pi_count] ) ||
			                ( matchType==2 && isGoodPion[pi_count] ) ||
			                ( matchType==3 && isGoodPion3d[pi_count] );

			if( !selected ) continue;

			if( bQ2_rec < 0 || bQ2_rec >= nQ2Bins ) continue;
			if( bxB_rec < 0 || bxB_rec >= nXbBins ) continue;

			double z_rec = pi.getZ();
			if( z_rec < .3 || z_rec >= 1. ) continue;

			// Fill N_rec
			TH3F* hR = (chargeIdx == 0) ? hRec_p : hRec_m;
			hR->Fill( xB_rec, Q2_rec, z_rec );

			// Fill N_match: rec pion has a matched gen pion in the same bin
			if( !isGoodGenPion[pi_count] ) continue;

			double theta_gen = pi_match[pi_count].get3Momentum().Theta() * rad_to_deg;
			double p_gen     = pi_match[pi_count].get3Momentum().Mag();

			if( theta_gen < theta_cut ) continue;

			bool gen_matching = true;
			if( matchType == 1 || matchType == 2 )
				gen_matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() );
			else if( matchType == 3 )
				gen_matching = anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 1 ) >= 0 &&
				               anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 2 ) >= 0;

			if( !gen_matching ) continue;
			if( anal.applyAcceptanceMap( e_MC->get3Momentum().Mag(), rad_to_deg*e_MC->get3Momentum().Phi(), rad_to_deg*e_MC->get3Momentum().Theta(), 0 ) < 0 ) continue;
			if( anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, chargeIdx + 1 ) < 0 ) continue;

			double z_gen = pi_match[pi_count].get3Momentum().Mag();   // use gen z
			z_gen = pi_match[pi_count].getZ();

			if( z_gen < .3 || z_gen >= 1. ) continue;
			if( bQ2_gen < 0 || bQ2_gen >= nQ2Bins ) continue;
			if( bxB_gen < 0 || bxB_gen >= nXbBins ) continue;

			// Check if rec bin == gen bin in all three dimensions
			int bZ_rec = hRec_p->GetZaxis()->FindBin( z_rec ) - 1;
			int bZ_gen = hRec_p->GetZaxis()->FindBin( z_gen ) - 1;

			if( bxB_rec == bxB_gen && bQ2_rec == bQ2_gen && bZ_rec == bZ_gen ){
				TH3F* hM = (chargeIdx == 0) ? hMatch_p : hMatch_m;
				hM->Fill( xB_rec, Q2_rec, z_rec );
			}
		}
	}

	// ── Generator tree loop ───────────────────────────────────────────────────
	// Fills hGen: all generated pions passing kinematic cuts.
	TTreeReader reader_gen(genChain);
	TTreeReaderValue<genElectron>  e_gen(  reader_gen, "e_gen");
	TTreeReaderArray<genPion>      pi_gen( reader_gen, "pi_gen");

	event_total = genChain->GetEntries();
	cout << "Processing gen tree (" << event_total << " events)...\n";

	while( reader_gen.Next() ){
		int event_count = reader_gen.GetCurrentEntry();
		if( event_count % 1000000 == 0 )
			cout << "  " << event_count << " / " << event_total << "\n";

		if( !anal.applyElectronKinematicCuts( *e_gen ) ) continue;

		double Q2 = e_gen->getQ2();
		double xB = e_gen->getXb();
		int bQ2 = (int)( (Q2 - Q2_min) / (Q2_max - Q2_min) * nQ2Bins );
		int bxB = (int)( (xB - xB_min) / (xB_max - xB_min) * nXbBins );
		if( bQ2 < 0 || bQ2 >= nQ2Bins ) continue;
		if( bxB < 0 || bxB >= nXbBins ) continue;

		int pi_count = -1;
		for( auto pi : pi_gen ){
			pi_count++;

			double phi   = pi.get3Momentum().Phi();
			double theta = pi.get3Momentum().Theta() * rad_to_deg;
			double p     = pi.get3Momentum().Mag();

			if( theta < theta_cut ) continue;
			if( !anal.applyPionKinematicCuts(pi) || pi.getMx() < 1.7 ) continue;

			int chargeIdx = (int)( pi.getCharge() < 1 );

			// Sector assignment (mirrors calcCorrections)
			double sector_i = -1;
			if( pi.getCharge() > 0 ){
				if     ( phi > -0.8  && phi < 0.25  ) sector_i = 1;
				else if( phi >= 0.25 && phi < 1.3   ) sector_i = 2;
				else if( phi >= 1.3  && phi <= 2.35 ) sector_i = 3;
				else if( phi > 2.35  || phi < -2.9  ) sector_i = 4;
				else if( phi > -2.9  && phi < -1.85 ) sector_i = 5;
				else                                   sector_i = 6;
			} else {
				if     ( phi > -0.25 && phi < 0.8   ) sector_i = 1;
				else if( phi >= 0.8  && phi < 1.85  ) sector_i = 2;
				else if( phi >= 1.85 && phi <= 2.9  ) sector_i = 3;
				else if( phi > 2.9   || phi < -2.4  ) sector_i = 4;
				else if( phi > -2.4  && phi < -1.25 ) sector_i = 5;
				else                                   sector_i = 6;
			}

			bool matching = true;
			if( matchType == 1 || matchType == 2 )
				matching = anal.acceptance_match_2d( theta, p, sector_i );
			else if( matchType == 3 )
				matching = anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 1 ) >= 0 &&
				           anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 2 ) >= 0;

			if( !matching ) continue;

			double z = pi.getZ();
			if( z < .3 || z >= 1. ) continue;

			TH3F* hG = (chargeIdx == 0) ? hGen_p : hGen_m;
			hG->Fill( xB, Q2, z );
		}
	}

	// ── Compute purity and stability ─────────────────────────────────────────
	cout << "Computing purity and stability...\n";

	for( int i = 1; i <= nXbBins; i++ ){
		for( int j = 1; j <= nQ2Bins; j++ ){
			for( int k = 1; k <= nZBins; k++ ){

				for( int c = 0; c < 2; c++ ){
					TH3F* hR = (c==0) ? hRec_p   : hRec_m;
					TH3F* hG = (c==0) ? hGen_p   : hGen_m;
					TH3F* hM = (c==0) ? hMatch_p  : hMatch_m;
					TH3F* hP = (c==0) ? hPurity_p : hPurity_m;
					TH3F* hS = (c==0) ? hStability_p : hStability_m;

					double n_rec   = hR->GetBinContent(i, j, k);
					double n_gen   = hG->GetBinContent(i, j, k);
					double n_match = hM->GetBinContent(i, j, k);

					if( n_rec > 0 ){
						hP->SetBinContent(i, j, k, n_match / n_rec);
						hP->SetBinError(  i, j, k, getUncertainty(n_match, n_rec));
					}
					if( n_gen > 0 ){
						hS->SetBinContent(i, j, k, n_match / n_gen);
						hS->SetBinError(  i, j, k, getUncertainty(n_match, n_gen));
					}
				}
			}
		}
	}

	// ── Write output ─────────────────────────────────────────────────────────
	outFile->cd();
	hRec_p->Write();      hRec_m->Write();
	hGen_p->Write();      hGen_m->Write();
	hMatch_p->Write();    hMatch_m->Write();
	hPurity_p->Write();   hPurity_m->Write();
	hStability_p->Write();hStability_m->Write();
	outFile->Close();

	cout << "Done. Output: " << outName << "\n";
	return 0;
}

double getUncertainty( double num, double den ){
	if( num <= 0 || den <= 0 ) return 0;
	return (num/den) * sqrt( 1./num + 1./den );
}
