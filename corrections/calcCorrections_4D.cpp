#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TNamed.h>
#include <TH1.h>
#include <TH3.h>
#include "electron.h"
#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
#include "constants.h"
#include <TChain.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TError.h"

using std::cerr;
using std::cout;

const int nXbBins = bins_xB;
const int nQ2Bins = bins_Q2;
const int nZBins = bins_Z;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );
double getVarVal(TString var,  electron  e, pion pi );
double getVarVal(TString var,  genElectron  e, genPion pi );


int main( int argc, char** argv){

	// Layout: [OutFile] [Matching] [Target] [VarName] [VarBins] [VarMin] [VarMax]
	//         [nRec] [nGen] [rec files...] [gen files...] (opt: matchingFile)
	if( argc < 10 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Out File] [Matching] [Target] [Var Name] [Var Bins] [Var Min] [Var Max] [nRec] [nGen] [rec files...] [gen files...] (opt: matching file)\n";
		cerr << "       Target: 0 = RGB/deuterium, 1 = RGA/proton\n";
		return -1;
	}

	TString outName  = argv[1];
	int matchType    = atoi(argv[2]);
	int rga          = atoi(argv[3]);
	TString var_name = argv[4];
	int bins_var     = atoi(argv[5]);
	double var_min   = atof(argv[6]);
	double var_max   = atof(argv[7]);
	int nRec         = atoi(argv[8]);
	int nGen         = atoi(argv[9]);
	const int firstRec = 10;
	const int firstGen = firstRec + nRec;
	TString matchingFile = "matchCut2D_10.root";
	if( argc > firstGen + nGen ) matchingFile = argv[firstGen + nGen];

	TChain * recChain = new TChain("ePi");
	cerr << "Rec files:\n";
	for( int i = 0; i < nRec; i++ ){
		recChain->Add( argv[firstRec + i] );
		cerr << "  " << argv[firstRec + i] << "\n";
	}

	TChain * genChain = new TChain("ePi");
	cerr << "Gen files:\n";
	for( int i = 0; i < nGen; i++ ){
		genChain->Add( argv[firstGen + i] );
		cerr << "  " << argv[firstGen + i] << "\n";
	}

	TFile * outFile = new TFile( outName, "RECREATE");

	TH1F * recHists[bins_var][nXbBins][nQ2Bins][2];
	TH1F * matchHists[bins_var][nXbBins][nQ2Bins][2];
	TH1F * genHists[bins_var][nXbBins][nQ2Bins][2];

	TH3F * binMigration_p[bins_var];
	TH3F * binMigration_m[bins_var];
	TH3F * binMigration_full[bins_var];

	TH3F * accCorrection_p[bins_var];
	TH3F * accCorrection_m[bins_var];
	TH3F * accCorrection_full[bins_var];

	TH3F * mcCorrection_p[bins_var];
	TH3F * mcCorrection_m[bins_var];
	TH3F * mcCorrection_full[bins_var];

	for(int pt = 0; pt < bins_var; pt++ ){
		binMigration_p[pt]    = new TH3F( Form("hBinMigrationP_%i", pt),  "hBinMigrationP",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		binMigration_m[pt]    = new TH3F( Form("hBinMigrationM_%i", pt),  "hBinMigrationM",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		binMigration_full[pt] = new TH3F( Form("hBinMigration_%i",  pt),  "hBinMigration",   nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);

		accCorrection_p[pt]    = new TH3F( Form("hAccCorrectionP_%i", pt), "hAccCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		accCorrection_m[pt]    = new TH3F( Form("hAccCorrectionM_%i", pt), "hAccCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		accCorrection_full[pt] = new TH3F( Form("hAccCorrection_%i",  pt), "hAccCorrection",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);

		mcCorrection_p[pt]    = new TH3F( Form("hMcCorrP_%i", pt), "hMcCorrP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		mcCorrection_m[pt]    = new TH3F( Form("hMcCorrM_%i", pt), "hMcCorrM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
		mcCorrection_full[pt] = new TH3F( Form("hMcCorr_%i",  pt), "hMcCorr",  nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1);
	}
	for(int pt = 0; pt < bins_var; pt++ ){
		for( int i = 0; i < nXbBins; i++ ){
			for( int j = 0; j < nQ2Bins; j++ ){

				recHists[pt][i][j][0]   = new TH1F( (TString)"recHist_P_"+Form("_%i_%i_%i",   pt, i, j), Form("recHist_P_%i_%i",   i, j), 14, .3, 1);
				matchHists[pt][i][j][0] = new TH1F( (TString)"matchHist_P_"+Form("_%i_%i_%i", pt, i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1);
				genHists[pt][i][j][0]   = new TH1F( (TString)"genHist_P_"+Form("_%i_%i_%i",   pt, i, j), Form("genHist_P_%i_%i",   i, j), 14, .3, 1);

				recHists[pt][i][j][1]   = new TH1F( (TString)"recHist_M_"+Form("_%i_%i_%i",   pt, i, j), Form("recHist_M_%i_%i",   i, j), 14, .3, 1);
				matchHists[pt][i][j][1] = new TH1F( (TString)"matchHist_M_"+Form("_%i_%i_%i", pt, i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1);
				genHists[pt][i][j][1]   = new TH1F( (TString)"genHist_M_"+Form("_%i_%i_%i",   pt, i, j), Form("genHist_M_%i_%i",   i, j), 14, .3, 1);

			}
		}
	}

	analyzer anal(0, -1);
	anal.setAnalyzerLevel(1);
	anal.setTarget( 0 );
	anal.loadCutValues(-1, 10.2);
	anal.loadSamplingFractionParams();
	anal.loadMatchingFunctions(matchingFile);
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");

	TTreeReader reader_rec(recChain);

	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");
	TTreeReaderArray<bool> isGoodGenPion(reader_rec, "isGoodGenPion");
	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<genElectron> e_MC(reader_rec, "e_gen");
	TTreeReaderArray<pion> pi_vec(reader_rec, "pi");
	TTreeReaderArray<genPion> pi_match(reader_rec, "pi_gen");
	//TTreeReaderValue<int> target(reader_rec, "target");

	int event_total = recChain->GetEntries();

	gErrorIgnoreLevel = kError;  // suppress TTreeReader/TChain tree-transition warnings
	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();
		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		//if( rga==1 && *target == 2112 ) continue;

		double Q2 = e->getQ2();
		double xB = e->getXb();
		double Q2_MC = e_MC->getQ2();
		double xB_MC = e_MC->getXb();

		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);

		int this_bin_Q2_MC = (int)( ( (Q2_MC - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB_MC = (int)( ( (xB_MC - xB_min)/(xB_max-xB_min) )*nXbBins);
		if( !anal.applyElectronDetectorCuts( *e ) ) continue;

		int pi_count = -1;
		for( auto pi : pi_vec ){
			pi_count++;
			int chargeIdx = (int)( pi.getCharge() < 1 );

			if( abs( pi.getPID() ) != 211 || abs(pi.getPID_eb()) == 2212 ){continue;}
			if( !anal.applyPionDetectorCuts( pi, *e ) ) continue;

			bool matching = true;
			double p_mom = pi.get3Momentum().Mag();

			if( matchType == 1 )
				matching = anal.applyAcceptanceMatching(pi, 2);
			else if( matchType == 3 )
				matching = anal.applyAcceptanceMatching(pi, 2) &&
				           anal.applyAcceptanceMap(p_mom, rad_to_deg*pi.get3Momentum().Phi(), rad_to_deg*pi.get3Momentum().Theta(), 1) >= 0 &&
				           anal.applyAcceptanceMap(p_mom, rad_to_deg*pi.get3Momentum().Phi(), rad_to_deg*pi.get3Momentum().Theta(), 2) >= 0;
			else if( matchType == 4 )
				matching = anal.applyAcceptanceMatching(pi, 2) &&
				           anal.applyAcceptanceMap(p_mom, rad_to_deg*pi.get3Momentum().Phi(), rad_to_deg*pi.get3Momentum().Theta(), 1, chargeIdx + 1) >= 0 &&
				           anal.applyAcceptanceMap(p_mom, rad_to_deg*pi.get3Momentum().Phi(), rad_to_deg*pi.get3Momentum().Theta(), 2, chargeIdx + 1) >= 0;

			if( matching ){
				if(anal.applyAcceptanceMap( e->get3Momentum().Mag(), rad_to_deg*e->get3Momentum().Phi(), rad_to_deg*e->get3Momentum().Theta(), 0 ) >= 0 &&
					anal.applyAcceptanceMap( p_mom, rad_to_deg*pi.get3Momentum().Phi(), rad_to_deg*pi.get3Momentum().Theta(), chargeIdx + 1 ) >= 0) {
					if( anal.applyElectronKinematicCuts( *e ) && anal.applyPionKinematicCuts( pi ) ){
						int this_bin_var = (int)( ( (getVarVal(var_name, *e, pi) - var_min)/(var_max - var_min) )*bins_var);
						if( this_bin_var >= 0 && this_bin_var < bins_var )
							recHists[this_bin_var][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
					}
				}
			}

			//Fill matched pions
			if( isGoodGenPion[pi_count] ){
				double theta_gen = pi_match[pi_count].get3Momentum().Theta()*rad_to_deg;
				double p_gen = pi_match[pi_count].get3Momentum().Mag();
				int chargeIdx_gen = (int)( pi_match[pi_count].getCharge() < 1 );

				double varval_MC = getVarVal(var_name, *e_MC, pi_match[pi_count]);
				int this_bin_var_MC = (int)( ( (varval_MC - var_min)/(var_max - var_min) )*bins_var);
				if( this_bin_var_MC < 0 || this_bin_var_MC >= bins_var ) continue;

				matching = true;
				if( matchType == 1 || matchType == 2 ){
					matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() );
				}
				else if( matchType == 3 ){
					matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() ) &&
					           anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 1 ) >= 0 &&
					           anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, 2 ) >= 0;
				}
				else if( matchType == 4 ){
					matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() ) &&
					           anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, chargeIdx_gen + 1 ) >= 0;
				}

				if( matching ){
					if(anal.applyAcceptanceMap( e_MC->get3Momentum().Mag(), rad_to_deg*e_MC->get3Momentum().Phi(), rad_to_deg*e_MC->get3Momentum().Theta(), 0 ) < 0) continue;
					if(anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(), theta_gen, chargeIdx_gen + 1 ) < 0 ) continue;
					if( !anal.applyElectronKinematicCuts( *e_MC ) ) continue;
					if( !anal.applyPionKinematicCuts( pi_match[pi_count] ) ) continue;
					matchHists[this_bin_var_MC][this_bin_xB_MC][this_bin_Q2_MC][chargeIdx_gen]->Fill( pi_match[pi_count].getZ() );
				}
			}
		}
	}

	TTreeReader reader_gen(genChain);
	TTreeReaderValue<genElectron> e_gen(reader_gen, "e_gen");
	TTreeReaderArray<genPion> pi_gen(reader_gen, "pi_gen");
	//TTreeReaderValue<int> target_gen(reader_gen, "target");

	event_total = genChain->GetEntries();

	gErrorIgnoreLevel = kError;
	while (reader_gen.Next()) {
		int event_count = reader_gen.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

		//if( rga == 1 && *target_gen != 2212 ) continue;

		// Apply electron kinematic cuts (W, Q2, xB, y, p_e, theta_e)
		// to match the cuts applied to reconstructed electrons via isGoodPion_no_acc
		if( !anal.applyElectronKinematicCuts( *e_gen ) ) continue;

		double Q2 = e_gen->getQ2();
		double xB = e_gen->getXb();

		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);

		int pi_count = -1;
		for( auto pi : pi_gen ){
			pi_count++;

			double phi = pi.get3Momentum().Phi();
			double theta = pi.get3Momentum().Theta()*rad_to_deg;
			double p = pi.get3Momentum().Mag();
			double charge = pi.getCharge();
			bool matching = true;
			int sector_i = -1;

			// Apply pion kinematic cuts (Mx, p_pi, z, theta_pi)
			// to match the cuts applied to reconstructed pions via isGoodPion_no_acc
			if( !anal.applyPionKinematicCuts(pi) || pi.getMx() < 1.7 ){continue;}

			int chargeIdx = (int)( pi.getCharge() < 1 );
			sector_i = -1;
			if( sector_i < 1 ){

				if( charge > 0 ){
					if( phi > -0.8 && phi < 0.25 ){ sector_i = 1; }
					else if( phi >= 0.25 && phi < 1.3 ){ sector_i = 2; }
					else if( phi >= 1.3 && phi <= 2.35 ){ sector_i = 3; }
					else if( phi > 2.35 || phi < -2.9  ){ sector_i = 4; }
					else if( phi > -2.9 && phi < -1.85){ sector_i = 5; }
					else{ sector_i = 6; }
				}
				if( charge < 0 ){
					if( phi > -0.25 && phi < 0.8 ){ sector_i = 1; }
					else if( phi >= 0.8 && phi < 1.85 ){ sector_i = 2; }
					else if( phi >= 1.85 && phi <= 2.9 ){ sector_i = 3; }
					else if( phi > 2.9 || phi < -2.4  ){ sector_i = 4; }
					else if( phi > -2.4 && phi < -1.25){ sector_i = 5; }
					else{ sector_i = 6; }
				}
			}

			if( matchType == 1 || matchType == 2 ){
				matching = anal.acceptance_match_2d( theta, p, sector_i );
			}
			else if( matchType == 3 ){
				matching = anal.acceptance_match_2d( theta, p, sector_i ) &&
				           anal.applyAcceptanceMap( p, rad_to_deg*pi.get3Momentum().Phi(), theta, 1 ) >= 0 &&
				           anal.applyAcceptanceMap( p, rad_to_deg*pi.get3Momentum().Phi(), theta, 2 ) >= 0;
			}
			else if( matchType == 4 ){
				matching = anal.acceptance_match_2d( theta, p, sector_i ) &&
				           anal.applyAcceptanceMap( p, rad_to_deg*pi.get3Momentum().Phi(), theta, chargeIdx + 1 ) >= 0;
			}
			else{ matching = true; }

			if( !matching ){ continue; }

			double varval = getVarVal(var_name, *e_gen, pi);
			int this_bin_var = (int)( ( (varval - var_min)/(var_max - var_min) )*bins_var);
			if( this_bin_var < 0 || this_bin_var >= bins_var ) continue;

			genHists[this_bin_var][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
		}
	}

	outFile->cd();

	for( int pt = 0; pt < bins_var; pt++ ){
		for( int i = 0; i < nXbBins; i++ ){
			for( int j = 0; j < nQ2Bins; j++ ){
				recHists[pt][i][j][0]->Write();
				recHists[pt][i][j][1]->Write();
				matchHists[pt][i][j][0]->Write();
				matchHists[pt][i][j][1]->Write();
				genHists[pt][i][j][0]->Write();
				genHists[pt][i][j][1]->Write();

				for( int k = 1; k <= nZBins; k++ ){

					double recBinPos   = recHists[pt][i][j][0]->GetBinContent(k);
					double recBinMin   = recHists[pt][i][j][1]->GetBinContent(k);
					double matchBinPos = matchHists[pt][i][j][0]->GetBinContent(k);
					double matchBinMin = matchHists[pt][i][j][1]->GetBinContent(k);
					double genBinPos   = genHists[pt][i][j][0]->GetBinContent(k);
					double genBinMin   = genHists[pt][i][j][1]->GetBinContent(k);

					double binMigrationPos = matchBinPos/recBinPos;
					double accCorrPos      = genBinPos/matchBinPos;
					double mcCorrPos       = genBinPos/recBinPos;
					double binMigrationMin = matchBinMin/recBinMin;
					double accCorrMin      = genBinMin/matchBinMin;
					double mcCorrMin       = genBinMin/recBinMin;

					double accCorr          = accCorrPos/accCorrMin;
					double binMigrationCorr = binMigrationPos/binMigrationMin;
					double mcCorr           = mcCorrPos/mcCorrMin;

					binMigration_p[pt]->SetBinContent( i+1, j+1, k, binMigrationPos );
					accCorrection_p[pt]->SetBinContent( i+1, j+1, k, accCorrPos );
					mcCorrection_p[pt]->SetBinContent( i+1, j+1, k, mcCorrPos );

					binMigration_m[pt]->SetBinContent( i+1, j+1, k, binMigrationMin );
					accCorrection_m[pt]->SetBinContent( i+1, j+1, k, accCorrMin );
					mcCorrection_m[pt]->SetBinContent( i+1, j+1, k, mcCorrMin );

					binMigration_full[pt]->SetBinContent( i+1, j+1, k, binMigrationCorr );
					accCorrection_full[pt]->SetBinContent( i+1, j+1, k, accCorr );
					mcCorrection_full[pt]->SetBinContent( i+1, j+1, k, mcCorr );

					binMigration_p[pt]->SetBinError( i+1, j+1, k, getUncertainty( matchBinPos, recBinPos ) );
					accCorrection_p[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinPos, matchBinPos ) );
					mcCorrection_p[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinPos, recBinPos ) );

					binMigration_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( matchBinMin, recBinMin ) );
					accCorrection_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, matchBinMin ) );
					mcCorrection_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, recBinMin ) );

					binMigration_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( binMigrationPos, binMigrationMin, getUncertainty( matchBinPos, recBinPos ), getUncertainty( matchBinMin, recBinMin ) ) );
					accCorrection_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( accCorrPos, accCorrMin, getUncertainty( genBinPos, matchBinPos ), getUncertainty( genBinMin, matchBinMin ) ) );
					mcCorrection_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( mcCorrPos, mcCorrMin, getUncertainty( genBinPos, recBinPos ), getUncertainty( genBinMin, recBinMin ) ) );
				}
			}
		}
		binMigration_p[pt]->Write();
		accCorrection_p[pt]->Write();
		mcCorrection_p[pt]->Write();

		binMigration_m[pt]->Write();
		accCorrection_m[pt]->Write();
		mcCorrection_m[pt]->Write();

		binMigration_full[pt]->Write();
		accCorrection_full[pt]->Write();
		mcCorrection_full[pt]->Write();
	}

	// Store binning metadata so consumers can self-configure
	TParameter<int>    nBinsPar("n4d_bins", bins_var);  nBinsPar.Write();
	TParameter<double> varMinPar("var_min",  var_min);   varMinPar.Write();
	TParameter<double> varMaxPar("var_max",  var_max);   varMaxPar.Write();
	TNamed varNameObj("var_name", var_name.Data());      varNameObj.Write();

	outFile->Close();

	return 1;
}

double getUncertainty( double num, double den ){
	return (num/den)*sqrt( pow( sqrt(num)/num , 2 ) + pow( sqrt(den)/den , 2 ) );
}

double getUncertainty( double num, double den, double num_unc, double den_unc ){
	return (num/den)*sqrt( pow( num_unc/num , 2 ) + pow( den_unc/den , 2 ) );
}


double getVarVal(TString var,  electron  e, pion pi ){
	if( var == "p_e" ) return e.get3Momentum().Mag();
	if( var == "theta_e" ) return e.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e.get3Momentum().Phi()*rad_to_deg;
	if( var == "W2" ) return e.getW2();
	if( var == "Q2" ) return e.getQ2();
	if( var == "xB" ) return e.getXb();
	if( var == "y" ) return e.getY();
	if( var == "sector" ) return e.getDC_sector();
	if( var == "sector_e" ) { double phi_deg = e.get3Momentum().Phi()*rad_to_deg + 20.; if (phi_deg < 0.) phi_deg += 360.; return (int)(phi_deg / 60.) + 1; }

	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) return pi.getPi_q().Phi()*rad_to_deg + 360*( (int) (pi.getPi_q().Phi() < 0 ) );
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return pi.getPi_q().Pt();
	if( var == "sector_pi" ) return pi.getDC_sector();
	if( var == "eta" ) return pi.getEta();

	return 0;
}


double getVarVal(TString var,  genElectron  e, genPion pi ){
	if( var == "p_e" ) return e.get3Momentum().Mag();
	if( var == "theta_e" ) return e.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e.get3Momentum().Phi()*rad_to_deg;
	if( var == "sector_e" ) { double phi_deg = e.get3Momentum().Phi()*rad_to_deg + 20.; if (phi_deg < 0.) phi_deg += 360.; return (int)(phi_deg / 60.) + 1; }
	if( var == "W2" ) return e.getW2();
	if( var == "Q2" ) return e.getQ2();
	if( var == "xB" ) return e.getXb();
	if( var == "y" ) return e.getY();

	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) return pi.getPi_q().Phi()*rad_to_deg + 360.*( (int) (pi.getPi_q().Phi() < 0 ) );
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return pi.getPi_q().Pt();
	if( var == "eta" ) return pi.getEta();

	return 0;
}
