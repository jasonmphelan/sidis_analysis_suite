#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TNamed.h>
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

const int nXbBins = bins_xB;
const int nQ2Bins = bins_Q2;
const int nZBins = bins_Z;

//const double xB_min = .1;
//const double xB_max = .6;

//const double Q2_min = 2;
//const double Q2_max = 8;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );
double getVarVal(TString var,  electron  e, pion pi );
double getVarVal(TString var,  genElectron  e, genPion pi );


int main( int argc, char** argv){
	
	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Rec File] [Gen File]  [Output File] [Matching]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2]  << "\n";
	
	TString inName_rec = argv[1];
	TString inName_gen = argv[2];
	TString outName = argv[3];
	int matchType = atoi( argv[4] );
	TString var_name = argv[5];
	int bins_var = atoi(argv[6]);
	double var_min = atof(argv[7]);
	double var_max = atof(argv[8]);


	TFile * outFile = new TFile( outName, "RECREATE");
	
	TFile * inFile_rec = new TFile( inName_rec );
	TTree * recChain = (TTree *)inFile_rec->Get("ePi");
	
	TFile * inFile_gen = new TFile( inName_gen );
	TTree * genChain = (TTree *)inFile_gen->Get("ePi");
	
	TH1F * recHists[bins_var][nXbBins][nQ2Bins][2];
	TH1F * matchHists[bins_var][nXbBins][nQ2Bins][2];
	TH1F * genHists[bins_var][nXbBins][nQ2Bins][2];
	
	TH3F * binMigration_p[bins_var]; 
	TH3F * binMigration_m[bins_var] ;
	TH3F * binMigration_full[bins_var] ;
	
	TH3F * accCorrection_p[bins_var] ;
	TH3F * accCorrection_m[bins_var];
	TH3F * accCorrection_full[bins_var]; 
	
	TH3F * mcCorrection_p[bins_var] ;
	TH3F * mcCorrection_m[bins_var];
	TH3F * mcCorrection_full[bins_var] ;

	for(int pt = 0; pt < bins_var; pt++ ){
		binMigration_p[pt] = new TH3F( Form("hBinMigrationP_%i", pt), "hBinMigrationP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		binMigration_m[pt] = new TH3F( Form("hBinMigrationM_%i", pt), "hBinMigrationM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		binMigration_full[pt] = new TH3F( Form("hBinMigration_%i", pt), "hBinMigration", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		
		accCorrection_p[pt] = new TH3F( Form("hAccCorrectionP_%i", pt), "hAccCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		accCorrection_m[pt] = new TH3F( Form("hAccCorrectionM_%i", pt), "hAccCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		accCorrection_full[pt] = new TH3F( Form("hAccCorrection_%i", pt), "hAccCorrection", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		
		mcCorrection_p[pt] = new TH3F( Form("hMcCorrP_%i", pt), "hKaonCorrP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		mcCorrection_m[pt] = new TH3F( Form("hMcCorrM_%i", pt), "hKaonCorrM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		mcCorrection_full[pt] = new TH3F( Form("hMcCorr_%i", pt), "hKaonCorr", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	}	
	for(int pt = 0; pt < bins_var; pt++ ){
		for( int i = 0; i < nXbBins; i++ ){
			for( int j = 0; j < nQ2Bins; j++ ){
				
				recHists[pt][i][j][0] = new TH1F( (TString)"recHist_P_"+Form("_%i_%i_%i", pt, i, j), Form("recHist_P_%i_%i", i, j), 14, .3, 1); 
				matchHists[pt][i][j][0] = new TH1F((TString)"matchHist_P_"+Form("_%i_%i_%i", pt, i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1); 
				genHists[pt][i][j][0] = new TH1F( (TString)"genHist_P_"+Form("_%i_%i_%i", pt, i, j), Form("genHist_P_%i_%i", i, j), 14, .3, 1); 
				
				recHists[pt][i][j][1] = new TH1F( (TString)"recHist_M_"+Form("_%i_%i_%i", pt, i, j), Form("recHist_M_%i_%i", i, j), 14, .3, 1); 
				matchHists[pt][i][j][1] = new TH1F( (TString)"matchHist_M_"+Form("_%i_%i_%i", pt, i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1); 
				genHists[pt][i][j][1] = new TH1F( (TString)"genHist_M_"+Form("_%i_%i_%i", pt, i, j), Form("genHist_M_%i_%i", i, j), 14, .3, 1); 
				
			}
		}
	}
	analyzer anal(0, -1);
	anal.setAnalyzerLevel(1);
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();	
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));
	
	TTreeReader reader_rec(recChain);
	
	
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");

	TTreeReaderArray<bool> isGoodGenPion(reader_rec, "isGoodGenPion");

	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<genElectron> e_MC(reader_rec, "e_gen");
	TTreeReaderArray<pion> pi_vec(reader_rec, "pi");
	TTreeReaderArray<genPion> pi_match(reader_rec, "pi_gen");
	
	//Define good event list and additional variables for output branches
	
	int event_total = recChain->GetEntries();
	
	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();
		
		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		
		double Q2 = e->getQ2();
		double xB = e->getXb();
		double Q2_MC = e_MC->getQ2();
		double xB_MC = e_MC->getXb();
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);
		
		int this_bin_Q2_MC = (int)( ( (Q2_MC - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB_MC = (int)( ( (xB_MC - xB_min)/(xB_max-xB_min) )*nXbBins);
		
		
		int pi_count = -1;
		for( auto pi : pi_vec ){
			pi_count++;
			int chargeIdx = (int)( pi.getCharge() < 1 );

			if( abs( pi.getPID() ) != 211 ){continue;}
	

			bool matching = true;
			if( (matchType==0 && isGoodPion_no_acc[pi_count])|| (matchType==2 && isGoodPion[pi_count]) || (matchType==3 && isGoodPion3d[pi_count]) ){//|| (matchType==3 && isGoodPion3d[pi_count])){
				double pt = pi.getPi_q().Pt();
				int this_bin_var = (int)( ( (getVarVal(var_name, *e, pi) - var_min)/(var_max - var_min) )*bins_var);
				//if( matchType == 2 ){ matching = !isGoodPion[pi_count]; }
				//else if( matchType == 3 ){ matching = !isGoodPion3d[pi_count]; }
				//else{ matching = false; }
			
				//if( matching ){ continue; }
			
				//Fill reco pions
				
				recHists[this_bin_var][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
			}
			
			//Fill matched pions
			if( isGoodGenPion[pi_count] ){
				double theta_gen = pi_match[pi_count].get3Momentum().Theta()*rad_to_deg;
				double p_gen = pi_match[pi_count].get3Momentum().Mag();
				int chargeIdx = (int)( pi_match[pi_count].getCharge() < 1 );
				double phi = pi_match[pi_count].get3Momentum().Phi();
				double theta = pi_match[pi_count].get3Momentum().Theta()*rad_to_deg;
				double p = pi_match[pi_count].get3Momentum().Mag();
				int sector_i = 0;
				
				if( var_name == "sector_pi"){
					//sector_i = anal.applyAcceptanceMap( p, phi*rad_to_deg, theta, chargeIdx + 1 ) + 1;
					if( sector_i < 1 ){
			
						if( chargeIdx == 0 ){    
							if( phi > -0.8 && phi < 0.25 ){ sector_i = 1; }//1.05
							else if( phi >= 0.25 && phi < 1.3 ){ sector_i = 2; }//1.05
							else if( phi >= 1.3 && phi <= 2.35 ){ sector_i = 3; }//1.05
							else if( phi > 2.35 || phi < -2.85  ){ sector_i = 4; }
							else if( phi > -2.85 && phi < -1.85){ sector_i = 5; }//1
							else if( phi > -1.85 && phi < .85){ sector_i = 6; }//1
					
						}
						if( chargeIdx == 1 ){
							if( phi > -0.2 && phi < 0.85 ){ sector_i = 1; } //1.05
							else if( phi >= 0.85 && phi < 1.9 ){ sector_i = 2; } //1.05
							else if( phi >= 1.9 && phi <= 2.95 ){ sector_i = 3; } //1.05
							else if( phi > 2.95 || phi < -2.25  ){ sector_i = 4; }
							else if( phi > -2.25 && phi < -1.25){ sector_i = 5; }
							else if( phi > -1.25 && phi < 1.45){ sector_i = 6; }//1
						}
					}	
				}
				double varval = (var_name == "sector_pi") ? sector_i : getVarVal(var_name, *e_MC, pi_match[pi_count]);
				int this_bin_var_MC = (int)( ( (varval - var_min)/(var_max - var_min) )*bins_var);
				if( this_bin_var_MC >= bins_var) continue;

				bool matching = true;
				if( matchType == 2 ){
					matching = anal.acceptance_match_2d( theta_gen, p_gen, pi.getDC_sector() );	
				}
				else if( matchType == 3 ){
					matching = (anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(),theta_gen, 1 ) >= 0 && 
									anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(),theta_gen, 2 ) >= 0 );

				}
				
				if( matching ){
					if(anal.applyAcceptanceMap( e_MC->get3Momentum().Mag(), rad_to_deg*e_MC->get3Momentum().Phi(), rad_to_deg*e_MC->get3Momentum().Theta(), 0 ) <0) continue;
					if(anal.applyAcceptanceMap( p_gen, rad_to_deg*pi_match[pi_count].get3Momentum().Phi(),theta_gen, chargeIdx + 1 ) < 0 ) continue;
			
				//if( pi_match[pi_count].getZ() < .3 || pi_match[pi_count].getZ() > 1 ){ continue; }
				//if( this_bin_Q2_MC < 0 || this_bin_Q2_MC >= nQ2Bins ){ continue; }
				//if( this_bin_xB_MC < 0 || this_bin_xB_MC >= nXbBins ){ continue; }
					matchHists[this_bin_var_MC][this_bin_xB_MC][this_bin_Q2_MC][chargeIdx]->Fill( pi_match[pi_count].getZ() );
				}
			}
		}
	}
	
	
	
	TTreeReader reader_gen(genChain);
	TTreeReaderValue<genElectron> e_gen(reader_gen, "e_gen");
	TTreeReaderArray<genPion> pi_gen(reader_gen, "pi_gen");
	
	
	event_total = genChain->GetEntries();
	
	while (reader_gen.Next()) {
		int event_count = reader_gen.GetCurrentEntry();
		
		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		
		double Q2 = e_gen->getQ2();
		double xB = e_gen->getXb();
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
		int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);
		
		//if(anal.applyAcceptanceMap( e_gen->get3Momentum().Mag(), rad_to_deg*e_gen->get3Momentum().Phi(), rad_to_deg*e_gen->get3Momentum().Theta(), 0 ) <0) continue;
		
		int pi_count = -1;
		for( auto pi : pi_gen ){
			pi_count++;
			
			double phi = pi.get3Momentum().Phi();
			double theta = pi.get3Momentum().Theta()*rad_to_deg;
			double p = pi.get3Momentum().Mag();
			double charge = pi.getCharge();
			bool matching = true;
			int sector_i = -1;	
			double pt = pi.getPi_q().Pt();

			

			if( !anal.applyPionKinematicCuts(pi) ){continue;}
			

			int chargeIdx = (int)( pi.getCharge() < 1 );
			//sector_i = anal.applyAcceptanceMap( p, phi*rad_to_deg, theta, chargeIdx + 1 ) + 1;
			if( sector_i < 1 ){
			
				if( chargeIdx == 0 ){    
					if( phi > -0.8 && phi < 0.25 ){ sector_i = 1; }//1.05
					else if( phi >= 0.25 && phi < 1.3 ){ sector_i = 2; }//1.05
					else if( phi >= 1.3 && phi <= 2.35 ){ sector_i = 3; }//1.05
					else if( phi > 2.35 || phi < -2.85  ){ sector_i = 4; }
					else if( phi > -2.85 && phi < -1.85){ sector_i = 5; }//1
					else if( phi > -1.85 && phi < .85){ sector_i = 6; }//1
			
				}
				if( chargeIdx == 1 ){
					if( phi > -0.2 && phi < 0.85 ){ sector_i = 1; } //1.05
					else if( phi >= 0.85 && phi < 1.9 ){ sector_i = 2; } //1.05
					else if( phi >= 1.9 && phi <= 2.95 ){ sector_i = 3; } //1.05
					else if( phi > 2.95 || phi < -2.25  ){ sector_i = 4; }
					else if( phi > -2.25 && phi < -1.25){ sector_i = 5; }
					else if( phi > -1.25 && phi < 1.45){ sector_i = 6; }//1
				}
			}	
		

			double varval = (var_name == "sector_pi") ? sector_i :  getVarVal(var_name, *e_gen, pi);
			int this_bin_var = (int)( ( (varval - var_min)/(var_max - var_min) )*bins_var);
			//if( matchType == 2 ){ matching = !isGoodPion_gen[pi_count]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d_MC[pi_count]; }
			//else{ matching = false; }
			
			if( matchType == 2 ){
				matching = anal.acceptance_match_2d( theta, p, sector_i );	
			}
			else if( matchType == 3 ){
				matching = (anal.applyAcceptanceMap( p, rad_to_deg*pi.get3Momentum().Phi(),theta, 1 ) >= 0 && 
								anal.applyAcceptanceMap( p, rad_to_deg*pi.get3Momentum().Phi(),theta, 2 ) >= 0 );

			}
			//else if( matchType == 3 ){
			//	phi = pi.get3Momentum().Phi()*rad_to_deg;
			//if( acceptance_match_3d( phi, theta, p, 0 ) && acceptance_match_3d( phi, theta, p, 1) ){
			//	if( anal.acceptance_match_3d_cont( phi, theta, p, 0 ) > -1
			//		&& anal.acceptance_match_3d_cont( phi, theta, p, 1 ) > -1   ){ matching = true; }
			//	else{ matching = false; }
			
			//}
			else{ matching = true ;}
			
			if( !matching ){ continue; }
			
			//Fill reco pions
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
					
					double recBinPos = recHists[pt][i][j][0]->GetBinContent(k);
					double recBinMin = recHists[pt][i][j][1]->GetBinContent(k);
					
					double matchBinPos = matchHists[pt][i][j][0]->GetBinContent(k);
					double matchBinMin = matchHists[pt][i][j][1]->GetBinContent(k);
					
					double genBinPos = genHists[pt][i][j][0]->GetBinContent(k);
					double genBinMin = genHists[pt][i][j][1]->GetBinContent(k);
					
					double binMigrationPos = matchBinPos/recBinPos;
					double accCorrPos = genBinPos/matchBinPos;				
					double mcCorrPos = genBinPos/recBinPos;
					
					double binMigrationMin = matchBinMin/recBinMin;
					double accCorrMin = genBinMin/matchBinMin;				
					double mcCorrMin = genBinMin/recBinMin;
					
					double accCorr = accCorrPos/accCorrMin;
					double binMigrationCorr = binMigrationPos/binMigrationMin;
					double mcCorr = mcCorrPos/mcCorrMin;
					
					
					
					
					binMigration_p[pt]->SetBinContent( i+1, j+1, k,  binMigrationPos );
					accCorrection_p[pt]->SetBinContent( i+1, j+1, k,  accCorrPos );
					mcCorrection_p[pt]->SetBinContent( i+1, j+1, k,  mcCorrPos );
					
					binMigration_m[pt]->SetBinContent( i+1, j+1, k, binMigrationMin );
					accCorrection_m[pt]->SetBinContent( i+1, j+1, k, accCorrMin );
					mcCorrection_m[pt]->SetBinContent( i+1, j+1, k, mcCorrMin );
					
					binMigration_full[pt]->SetBinContent( i+1, j+1, k, binMigrationCorr );
					accCorrection_full[pt]->SetBinContent( i+1, j+1, k, accCorr );
					mcCorrection_full[pt]->SetBinContent( i+1, j+1, k, mcCorr );
					
					binMigration_p[pt]->SetBinError( i+1, j+1, k,  getUncertainty( matchBinPos, recBinPos ) );
					accCorrection_p[pt]->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, matchBinPos ) );
					mcCorrection_p[pt]->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, recBinPos ) );
					
					binMigration_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( matchBinMin, recBinMin ) );
					accCorrection_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, matchBinMin) );
					mcCorrection_m[pt]->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, recBinMin) );
					
					binMigration_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( binMigrationPos, binMigrationMin, getUncertainty( matchBinPos, recBinPos ), getUncertainty( matchBinMin, recBinMin  ) ) );
					accCorrection_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( accCorrPos, accCorrMin, getUncertainty(genBinPos, matchBinPos ),getUncertainty(genBinMin, matchBinMin )  ) );
					mcCorrection_full[pt]->SetBinError( i+1, j+1, k, getUncertainty( mcCorrPos, mcCorrMin, getUncertainty(genBinPos, recBinPos ),getUncertainty(genBinMin, recBinMin )  ) );
					
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
	TParameter<int>   nBinsPar("n4d_bins", bins_var);   nBinsPar.Write();
	TParameter<double> varMinPar("var_min",  var_min);   varMinPar.Write();
	TParameter<double> varMaxPar("var_max",  var_max);   varMaxPar.Write();
	TNamed varNameObj("var_name", var_name.Data());      varNameObj.Write();

	outFile->Close();
	//file_1->Close();
	//file_2->Close();
	
	return 1;
}

double getUncertainty( double num, double den ){
	return (num/den)*sqrt( pow( sqrt(num)/num , 2 ) + pow( sqrt(den)/den , 2 ) );
	//return sqrt( pow( sqrt(num)/num , 2 ) + pow( sqrt(den)/den , 2 ) );
	
}

double getUncertainty( double num, double den, double num_unc, double den_unc ){
	return (num/den)*sqrt( pow( num_unc/num , 2 ) + pow( den_unc/den , 2 ) );
	//return sqrt( pow( num_unc/num , 2 ) + pow( den_unc/den , 2 ) );
	
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
	if( var == "phi_q" ) {
		return pi.getPi_q().Phi()*rad_to_deg + 360*( (int) (pi.getPi_q().Phi() < 0 ) ) ;
	}
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return  pi.getPi_q().Pt();
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
	if( var == "phi_q" ) {
		return pi.getPi_q().Phi()*rad_to_deg + 360.*( (int) (pi.getPi_q().Phi() < 0 ) ) ;
	}
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return  pi.getPi_q().Pt();
	if( var == "eta" ) return pi.getEta();

	return 0;
}


