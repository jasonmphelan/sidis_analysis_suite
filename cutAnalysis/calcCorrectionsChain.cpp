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
#include <TH3.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "electron.h"
#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
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

const int nXbBins = bins_xB;
const int nQ2Bins = bins_Q2;
const int nZBins  = bins_Z;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );

// Look up the mcCorrection value at (xB, Q2, z) from the in-memory TH3F.
// Returns 1 if the bin is empty or the value is non-finite.
double lookupCorr( TH3F* h, double xB, double Q2, double z ){
	if( !h ) return 1.0;
	int bx = h->GetXaxis()->FindBin(xB);
	int bq = h->GetYaxis()->FindBin(Q2);
	int bz = h->GetZaxis()->FindBin(z);
	double val = h->GetBinContent(bx, bq, bz);
	return (val > 0 && std::isfinite(val)) ? val : 1.0;
}

int main( int argc, char** argv){

	if( argc < 9 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# Files] [Beam Energy] [Run Type] [Output File] [Matching]\n";
		cerr << "       [Rec Skim Base] [Gen Skim Base] [Data Skim Base] [Target]\n";
		cerr << "  # Files        : 0 = all files, N = cap at N\n";
		cerr << "  Run Type       : 1=MC detector, 2=MC generator (for rec/gen chains)\n";
		cerr << "  Output File    : base name (no extension); writes _pip.txt and _pim.txt\n";
		cerr << "  Matching       : 0=no acc, 2=2D match, 3=3D match\n";
		cerr << "  Data Skim Base : base path for data skims (RunType=0 used automatically)\n";
		cerr << "  Target         : 0 = RGB/deuterium, 1 = RGA/proton\n";
		return -1;
	}

	int     nFiles    = atoi(argv[1]);
	double  Ebeam     = atof(argv[2]);
	TString outName   = argv[3];
	int     matchType = atoi(argv[4]);
	TString recBase   = argv[5];
	TString genBase   = argv[6];
	TString dataBase  = argv[7];
	int     target    = atoi(argv[8]); // 0 = RGB/deuterium, 1 = RGA/proton

	cerr << "Rec base  : " << recBase  << "\n";
	cerr << "Gen base  : " << genBase  << "\n";
	cerr << "Data base : " << dataBase << "\n";

	////////////////////////////////////////////////////
	// Build MC rec and gen TChains via reader
	////////////////////////////////////////////////////
	TChain * recChain  = new TChain("ePi");
	TChain * genChain  = new TChain("ePi");
	TChain * dataChain = new TChain("ePi");

	recChain->Add(recBase);
	genChain->Add(genBase);
	dataChain->Add(dataBase);

	cout << "Rec  chain entries : " << recChain->GetEntries()  << "\n";
	cout << "Gen  chain entries : " << genChain->GetEntries()  << "\n";
	cout << "Data chain entries : " << dataChain->GetEntries() << "\n";

	////////////////////////////////////////////////////
	// Allocate intermediate and correction histograms
	////////////////////////////////////////////////////

	// Intermediate z-distributions: [xB][Q2][charge]
	TH1F * recHists[nXbBins][nQ2Bins][2];
	TH1F * matchHists[nXbBins][nQ2Bins][2];
	TH1F * genHists[nXbBins][nQ2Bins][2];

	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			recHists[i][j][0]   = new TH1F( Form("recHist_P_%i_%i",   i, j), Form("recHist_P_%i_%i",   i, j), 14, .3, 1 );
			matchHists[i][j][0] = new TH1F( Form("matchHist_P_%i_%i", i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1 );
			genHists[i][j][0]   = new TH1F( Form("genHist_P_%i_%i",   i, j), Form("genHist_P_%i_%i",   i, j), 14, .3, 1 );

			recHists[i][j][1]   = new TH1F( Form("recHist_M_%i_%i",   i, j), Form("recHist_M_%i_%i",   i, j), 14, .3, 1 );
			matchHists[i][j][1] = new TH1F( Form("matchHist_M_%i_%i", i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1 );
			genHists[i][j][1]   = new TH1F( Form("genHist_M_%i_%i",   i, j), Form("genHist_M_%i_%i",   i, j), 14, .3, 1 );
		}
	}

	// In-memory correction TH3Fs: axes are (xB, Q2, z)
	TH3F * mcCorrection_p = new TH3F( "hMcCorrectionP", "hMcCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );
	TH3F * mcCorrection_m = new TH3F( "hMcCorrectionM", "hMcCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );

	// (Keep full set for optional diagnostics / future use)
	TH3F * binMigration_p    = new TH3F( "hBinMigrationP",  "hBinMigrationP",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );
	TH3F * binMigration_m    = new TH3F( "hBinMigrationM",  "hBinMigrationM",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );
	TH3F * accCorrection_p   = new TH3F( "hAccCorrectionP", "hAccCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );
	TH3F * accCorrection_m   = new TH3F( "hAccCorrectionM", "hAccCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1 );

	////////////////////////////////////////////////////
	// Load acceptance map and matching functions
	////////////////////////////////////////////////////
	analyzer anal(0, -1);
	anal.setAnalyzerLevel(1);
	anal.setTarget( target );
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root" );
	anal.loadSamplingFractionParams();
	////////////////////////////////////////////////////
	// MC rec loop: fill recHists and matchHists
	////////////////////////////////////////////////////
	cout << "--- MC rec loop ---\n";
	TTreeReader reader_rec(recChain);


	TTreeReaderValue<electron>    e_rec(reader_rec,      "e");
	TTreeReaderValue<genElectron> e_MC(reader_rec,       "e_gen");
	TTreeReaderArray<pion>        pi_vec(reader_rec,     "pi");
	TTreeReaderArray<genPion>     pi_match(reader_rec,   "pi_gen");

	int event_total = recChain->GetEntries();

	while( reader_rec.Next() ){
		int event_count = reader_rec.GetCurrentEntry();
		if( event_count%100000 == 0 ){
			cout << "MC rec: " << event_count << " / " << event_total << "\n";
		}

		double Q2    = e_rec->getQ2();
		double xB    = e_rec->getXb();
		double Q2_MC = e_MC->getQ2();
		double xB_MC = e_MC->getXb();

		int bin_Q2    = (int)( ( (Q2    - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins );
		int bin_xB    = (int)( ( (xB    - xB_min)/(xB_max-xB_min) )*nXbBins );
		int bin_Q2_MC = (int)( ( (Q2_MC - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins );
		int bin_xB_MC = (int)( ( (xB_MC - xB_min)/(xB_max-xB_min) )*nXbBins );

		int pi_count = -1;
		for( auto pi : pi_vec ){
			pi_count++;
			if( abs( pi.getPID() ) != 211 ){ continue; }

			int chargeIdx = (int)( pi.getCharge() < 1 );

			// Reconstructed pions

			if( anal.applyElectronDetectorCuts( *e_rec ) && anal.applyElectronKinematicCuts( *e_rec ) && anal.applyPionDetectorCuts( pi, *e_rec ) && anal.applyPionKinematicCuts( pi ) ){
				double p_pi = pi.get3Momentum().Mag();
				if( anal.applyAcceptanceMap( e_rec->get3Momentum().Mag(),
				                             rad_to_deg*e_rec->get3Momentum().Phi(),
				                             rad_to_deg*e_rec->get3Momentum().Theta(), 0 ) >= 0 &&
				    anal.applyAcceptanceMap( p_pi,
				                             rad_to_deg*pi.get3Momentum().Phi(),
				                             rad_to_deg*pi.get3Momentum().Theta(), chargeIdx + 1 ) >= 0 ){
					if( bin_xB >= 0 && bin_xB < nXbBins && bin_Q2 >= 0 && bin_Q2 < nQ2Bins )
						recHists[bin_xB][bin_Q2][chargeIdx]->Fill( pi.getZ() );
				}
			}

			// MC-truth-matched pions
			if( anal.applyElectronKinematicCuts( *e_MC ) && anal.applyPionKinematicCuts( pi_match[pi_count] ) ){
				double theta_gen = pi_match[pi_count].get3Momentum().Theta()*rad_to_deg;
				double p_gen     = pi_match[pi_count].get3Momentum().Mag();

				bool matching = true;
				if( matchType == 2 ){
					matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() );
				}
				else if( matchType == 3 ){
					matching = ( anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 1 ) >= 0 &&
					             anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 2 ) >= 0 );
				}

				if( matching ){
					if( anal.applyAcceptanceMap( e_MC->get3Momentum().Mag(),
					                             rad_to_deg*e_MC->get3Momentum().Phi(),
					                             rad_to_deg*e_MC->get3Momentum().Theta(), 0 ) < 0 ){ continue; }
					if( anal.applyAcceptanceMap( p_gen,
					                             rad_to_deg*pi_match[pi_count].get3Momentum().Phi(),
					                             theta_gen, chargeIdx + 1 ) < 0 ){ continue; }

					if( bin_xB_MC >= 0 && bin_xB_MC < nXbBins && bin_Q2_MC >= 0 && bin_Q2_MC < nQ2Bins )
						matchHists[bin_xB_MC][bin_Q2_MC][chargeIdx]->Fill( pi_match[pi_count].getZ() );
				}
			}
		}
	}

	////////////////////////////////////////////////////
	// MC gen loop: fill genHists
	////////////////////////////////////////////////////
	cout << "--- MC gen loop ---\n";
	TTreeReader reader_gen(genChain);
	TTreeReaderValue<genElectron> e_gen(reader_gen,  "e_gen");
	TTreeReaderArray<genPion>     pi_gen(reader_gen, "pi_gen");

	event_total = genChain->GetEntries();

	while( reader_gen.Next() ){
		int event_count = reader_gen.GetCurrentEntry();
		if( event_count%1000000 == 0 ){
			cout << "MC gen: " << event_count << " / " << event_total << "\n";
		}

		double Q2 = e_gen->getQ2();
		double xB = e_gen->getXb();
		int bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins );
		int bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins );

		int pi_count = -1;
		for( auto pi : pi_gen ){
			pi_count++;

			double phi      = pi.get3Momentum().Phi();
			double theta    = pi.get3Momentum().Theta()*rad_to_deg;
			double p        = pi.get3Momentum().Mag();
			int    chargeIdx = (int)( pi.getCharge() < 1 );

			bool   matching  = true;
			double sector_i  = -1;//anal.applyAcceptanceMap( p, phi*rad_to_deg, theta, chargeIdx + 1 ) + 1;
			if( !anal.applyElectronKinematicCuts( *e_gen ) || !anal.applyPionKinematicCuts( pi ) ) continue;

			if( sector_i < 1 ){
				if( pi.getCharge() > 0 ){    
					if( phi > -0.8 && phi < 0.25 ){ sector_i = 1; }
					else if( phi >= 0.25 && phi < 1.3 ){ sector_i = 2; }
					else if( phi >= 1.3 && phi <= 2.35 ){ sector_i = 3; }
					else if( phi > 2.35 || phi < -2.9  ){ sector_i = 4; }
					else if( phi > -2.9 && phi < -1.85){ sector_i = 5; }
					else{ sector_i = 6; }
				}
				if( pi.getCharge() < 0 ){
					if( phi > -0.25 && phi < 0.8 ){ sector_i = 1; }
					else if( phi >= 0.8 && phi < 1.85 ){ sector_i = 2; }
					else if( phi >= 1.85 && phi <= 2.9 ){ sector_i = 3; }
					else if( phi > 2.9 || phi < -2.4  ){ sector_i = 4; }
					else if( phi > -2.4 && phi < -1.25){ sector_i = 5; }
					else{ sector_i = 6; }
				}
			}

			if( matchType == 2 ){
				matching = anal.acceptance_match_2d( theta, p, sector_i );
			}
			else if( matchType == 3 ){
				matching = ( anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 1 ) >= 0 &&
				             anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 2 ) >= 0 );
			}
			else{ matching = true; }

			if( !matching ){ continue; }

			if( bin_xB >= 0 && bin_xB < nXbBins && bin_Q2 >= 0 && bin_Q2 < nQ2Bins )
				genHists[bin_xB][bin_Q2][chargeIdx]->Fill( pi.getZ() );
		}
	}

	////////////////////////////////////////////////////
	// Build in-memory correction TH3Fs
	////////////////////////////////////////////////////
	cout << "--- Building correction maps ---\n";

	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			for( int k = 1; k <= nZBins; k++ ){

				double recPos   = recHists[i][j][0]->GetBinContent(k);
				double recMin   = recHists[i][j][1]->GetBinContent(k);
				double matchPos = matchHists[i][j][0]->GetBinContent(k);
				double matchMin = matchHists[i][j][1]->GetBinContent(k);
				double genPos   = genHists[i][j][0]->GetBinContent(k);
				double genMin   = genHists[i][j][1]->GetBinContent(k);

				double mcCorrPos = (recPos  > 0) ? genPos/recPos   : 0;
				double mcCorrMin = (recMin  > 0) ? genMin/recMin   : 0;
				double accCorrPos = (matchPos > 0) ? genPos/matchPos : 0;
				double accCorrMin = (matchMin > 0) ? genMin/matchMin : 0;
				double binMigPos  = (recPos  > 0) ? matchPos/recPos : 0;
				double binMigMin  = (recMin  > 0) ? matchMin/recMin : 0;

				mcCorrection_p->SetBinContent(  i+1, j+1, k, mcCorrPos );
				mcCorrection_m->SetBinContent(  i+1, j+1, k, mcCorrMin );
				accCorrection_p->SetBinContent( i+1, j+1, k, accCorrPos );
				accCorrection_m->SetBinContent( i+1, j+1, k, accCorrMin );
				binMigration_p->SetBinContent(  i+1, j+1, k, binMigPos );
				binMigration_m->SetBinContent(  i+1, j+1, k, binMigMin );

				mcCorrection_p->SetBinError(  i+1, j+1, k, getUncertainty(genPos,  recPos)   );
				mcCorrection_m->SetBinError(  i+1, j+1, k, getUncertainty(genMin,  recMin)   );
				accCorrection_p->SetBinError( i+1, j+1, k, getUncertainty(genPos,  matchPos) );
				accCorrection_m->SetBinError( i+1, j+1, k, getUncertainty(genMin,  matchMin) );
				binMigration_p->SetBinError(  i+1, j+1, k, getUncertainty(matchPos, recPos)  );
				binMigration_m->SetBinError(  i+1, j+1, k, getUncertainty(matchMin, recMin)  );
			}
		}
	}

	////////////////////////////////////////////////////
	// Data loop: apply in-memory corrections
	////////////////////////////////////////////////////
	cout << "--- Data loop ---\n";

	// Output ratio histograms: [Q2][xB][charge]
	TH1F * hZ[nQ2Bins][nXbBins][2];
	TString charge_str[2] = {"", "Pim"};
	for( int i = 0; i < nQ2Bins; i++ ){
		for( int j = 0; j < nXbBins; j++ ){
			for( int c = 0; c < 2; c++ ){
				hZ[i][j][c] = new TH1F( "hRatio" + charge_str[c] + Form("_%i_%i", i+1, j+1),
				                         "hRatio" + charge_str[c] + Form("_%i_%i", i+1, j+1),
				                         nZBins, .3, 1 );
			}
		}
	}

	TTreeReader reader_data(dataChain);
	TTreeReaderValue<electron>    e_data(reader_data, "e");
	TTreeReaderArray<pion>        pi_data(reader_data, "pi");
	//TTreeReaderArray<bool>        isGoodPion_data(reader_data,        "isGoodPion");
	//TTreeReaderArray<bool>        isGoodPion_no_acc_data(reader_data, "isGoodPion_no_acc");
	//TTreeReaderArray<bool>        isGoodPion3d_data(reader_data,      "isGoodPion_3d");

	event_total = dataChain->GetEntries();

	while( reader_data.Next() ){
		int event_count = reader_data.GetCurrentEntry();
		if( event_count%100000 == 0 ){
			cout << "Data: " << event_count << " / " << event_total << "\n";
		}

		double xB = e_data->getXb();
		double Q2 = e_data->getQ2();

		// Acceptance map on electron
		TVector3 e_mom = e_data->get3Momentum();
		if( anal.applyAcceptanceMap( e_mom.Mag(), rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) < 0 ){ continue; }

		int bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins );
		int bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins );

		if( bin_Q2 < 0 || bin_Q2 >= nQ2Bins ){ continue; }
		if( bin_xB < 0 || bin_xB >= nXbBins ){ continue; }


		if( !anal.applyElectronDetectorCuts( *e_data ) || !anal.applyElectronKinematicCuts( *e_data ) ) continue;
					

		for( int i = 0; i < (int)(pi_data.end() - pi_data.begin()); i++ ){
			cout<<"enter pion loop \n";
			int chargeIdx = (int)( pi_data[i].getCharge() < 1 );
			double p_pi   = pi_data[i].get3Momentum().Mag();
			double z      = pi_data[i].getZ();
	
			// Additional matching cut
			bool matching = true;
			if( matchType == 2 ){ matching = anal.applyAcceptanceMatching(pi_data[i], 2); }
			else if( matchType == 3 ){
				matching = anal.applyAcceptanceMap( p_pi, rad_to_deg*pi_data[i].get3Momentum().Phi(),
													rad_to_deg*pi_data[i].get3Momentum().Theta(), 1 ) >= 0 &&
						anal.applyAcceptanceMap( p_pi, rad_to_deg*pi_data[i].get3Momentum().Phi(),
													rad_to_deg*pi_data[i].get3Momentum().Theta(), 2 ) >= 0;
			}
			if( !matching ){ continue; }

			if( !anal.applyPionDetectorCuts( pi_data[i], *e_data ) || !anal.applyPionKinematicCuts( pi_data[i] ) ) continue;

			// Acceptance map on pion
			if( anal.applyAcceptanceMap( p_pi, rad_to_deg*pi_data[i].get3Momentum().Phi(),
			                              rad_to_deg*pi_data[i].get3Momentum().Theta(), chargeIdx + 1 ) < 0 ){ continue; }

			

			// Apply MC correction weight from in-memory TH3F
			TH3F * corrHist = (chargeIdx == 0) ? mcCorrection_p : mcCorrection_m;
			double weight   = lookupCorr( corrHist, xB, Q2, z );
			hZ[bin_Q2][bin_xB][chargeIdx]->Fill( z, weight );
		}
	}

	////////////////////////////////////////////////////
	// Write output text files (same format as cutSensitivity.cpp)
	////////////////////////////////////////////////////
	cout << "--- Writing output ---\n";

	std::ofstream txtFile_pip, txtFile_pim;
	txtFile_pip.open( outName + "_pip.txt" );
	txtFile_pim.open( outName + "_pim.txt" );

	// Header: xB  Q2  1  2  ...  nZBins
	txtFile_pip << "xB\tQ2\t";
	txtFile_pim << "xB\tQ2\t";
	for( int z = 1; z <= nZBins; z++ ){
		txtFile_pip << z << "\t";
		txtFile_pim << z << "\t";
	}
	txtFile_pip << "\n";
	txtFile_pim << "\n";

	// One row per (xB, Q2) bin — outer loop xB, inner Q2 (matches cutSensitivity.cpp)
	for( int j = 1; j <= nXbBins; j++ ){
		for( int i = 1; i <= nQ2Bins; i++ ){

			txtFile_pip << j << "\t" << i;
			txtFile_pim << j << "\t" << i;

			for( int z = 1; z <= nZBins; z++ ){
				txtFile_pip << "\t" << hZ[i-1][j-1][0]->GetBinContent(z);
				txtFile_pim << "\t" << hZ[i-1][j-1][1]->GetBinContent(z);
			}
			txtFile_pip << "\n";
			txtFile_pim << "\n";
		}
	}

	txtFile_pip.close();
	txtFile_pim.close();
	cout << "Done.\n";

	return 1;
}

double getUncertainty( double num, double den ){
	if( den == 0 || num == 0 ) return 0;
	return (num/den)*sqrt( pow( sqrt(num)/num, 2 ) + pow( sqrt(den)/den, 2 ) );
}

double getUncertainty( double num, double den, double num_unc, double den_unc ){
	if( den == 0 || num == 0 ) return 0;
	return (num/den)*sqrt( pow( num_unc/num, 2 ) + pow( den_unc/den, 2 ) );
}
