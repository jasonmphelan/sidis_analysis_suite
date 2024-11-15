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

const int nXbBins = bins_xB;
const int nQ2Bins = bins_Q2;
const int nZBins = bins_Z;

//const double xB_min = .1;
//const double xB_max = .6;

//const double Q2_min = 2;
//const double Q2_max = 8;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );

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

	TFile * outFile = new TFile( outName, "RECREATE");

	TFile * inFile_rec = new TFile( inName_rec );
	TTree * recChain = (TTree *)inFile_rec->Get("ePi");
	
	TFile * inFile_gen = new TFile( inName_gen );
	TTree * genChain = (TTree *)inFile_gen->Get("ePi");

	TH1F * recHists[nXbBins][nQ2Bins][2];
	TH1F * matchHists[nXbBins][nQ2Bins][2];
	TH1F * genHists[nXbBins][nQ2Bins][2];

	TH3F * binMigration_p = new TH3F( "hBinMigrationP", "hBinMigrationP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * binMigration_m = new TH3F( "hBinMigrationM", "hBinMigrationM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * binMigration_full = new TH3F( "hBinMigration", "hBinMigration", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	
	TH3F * accCorrection_p = new TH3F( "hAccCorrectionP", "hAccCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * accCorrection_m = new TH3F( "hAccCorrectionM", "hAccCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * accCorrection_full = new TH3F( "hAccCorrection", "hAccCorrection", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	
	TH3F * mcCorrection_p = new TH3F( "hMcCorrectionP", "hMcCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * mcCorrection_m = new TH3F( "hMcCorrectionM", "hMcCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * mcCorrection_full = new TH3F( "hMcCorrection", "hMcCorrection", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	

	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
		
			recHists[i][j][0] = new TH1F( (TString)"recHist_P_"+Form("_%i_%i", i, j), Form("recHist_P_%i_%i", i, j), 14, .3, 1); 
			matchHists[i][j][0] = new TH1F((TString)"matchHist_P_"+Form("_%i_%i", i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1); 
			genHists[i][j][0] = new TH1F( (TString)"genHist_P_"+Form("_%i_%i", i, j), Form("genHist_P_%i_%i", i, j), 14, .3, 1); 
			
			recHists[i][j][1] = new TH1F( (TString)"recHist_M_"+Form("_%i_%i", i, j), Form("recHist_M_%i_%i", i, j), 14, .3, 1); 
			matchHists[i][j][1] = new TH1F( (TString)"matchHist_M_"+Form("_%i_%i", i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1); 
			genHists[i][j][1] = new TH1F( (TString)"genHist_M_"+Form("_%i_%i", i, j), Form("genHist_M_%i_%i", i, j), 14, .3, 1); 
			
		}
	}

        TTreeReader reader_rec(recChain);

        
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
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
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[pi_count]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[pi_count]; }
			else{ matching = false; }

			if( matching ){ continue; }

			int chargeIdx = (int)( pi.getCharge() < 1 );
		
			//Fill reco pions
			recHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
			

			//Fill matched pions
			if( pi_match[pi_count].getZ() < .3 || pi_match[pi_count].getZ() > 1 ){ continue; }
			if( this_bin_Q2_MC < 0 || this_bin_Q2_MC >= nQ2Bins ){ continue; }
			if( this_bin_xB_MC < 0 || this_bin_xB_MC >= nXbBins ){ continue; }
			matchHists[this_bin_xB_MC][this_bin_Q2_MC][chargeIdx]->Fill( pi_match[pi_count].getZ() );
			
		}
	}
				

	analyzer anal(0, -1);
	anal.setAnalyzerLevel(0);
	anal.loadCutValues(-1, 10.2);

	
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


		int pi_count = -1;
		for( auto pi : pi_gen ){
			pi_count++;

			double phi = pi.get3Momentum().Phi();
			double theta = pi.get3Momentum().Theta();
			double p = pi.get3Momentum().Mag();
			double charge = pi.getCharge();
			bool matching = true;
			double sector_i = -1;	
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

			//if( matchType == 2 ){ matching = !isGoodPion_gen[pi_count]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d_MC[pi_count]; }
			//else{ matching = false; }

			if( !anal.acceptance_match_2d( theta*rad_to_deg, p, sector_i )  ){ continue; }
			int chargeIdx = (int)( pi.getCharge() < 1 );
			
			//Fill reco pions
			genHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );

		}
	}

	

	outFile->cd();


	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			recHists[i][j][0]->Write();
			recHists[i][j][1]->Write();
			matchHists[i][j][0]->Write();
			matchHists[i][j][1]->Write();
			genHists[i][j][0]->Write();
			genHists[i][j][1]->Write();
			
			for( int k = 1; k <= nZBins; k++ ){

				double recBinPos = recHists[i][j][0]->GetBinContent(k);
				double recBinMin = recHists[i][j][1]->GetBinContent(k);

				double matchBinPos = matchHists[i][j][0]->GetBinContent(k);
				double matchBinMin = matchHists[i][j][1]->GetBinContent(k);

				double genBinPos = genHists[i][j][0]->GetBinContent(k);
				double genBinMin = genHists[i][j][1]->GetBinContent(k);
	
				double binMigrationPos = matchBinPos/recBinPos;
				double accCorrPos = genBinPos/matchBinPos;				
				double mcCorrPos = genBinPos/recBinPos;

				double binMigrationMin = matchBinMin/recBinMin;
				double accCorrMin = genBinMin/matchBinMin;				
				double mcCorrMin = genBinMin/recBinMin;
	
				double accCorr = accCorrPos/accCorrMin;
				double binMigrationCorr = binMigrationPos/binMigrationMin;
				double mcCorr = mcCorrPos/mcCorrMin;
				
				


				binMigration_p->SetBinContent( i+1, j+1, k,  binMigrationPos );
				accCorrection_p->SetBinContent( i+1, j+1, k,  accCorrPos );
				mcCorrection_p->SetBinContent( i+1, j+1, k,  mcCorrPos );
			
				binMigration_m->SetBinContent( i+1, j+1, k, binMigrationMin );
				accCorrection_m->SetBinContent( i+1, j+1, k, accCorrMin );
				mcCorrection_m->SetBinContent( i+1, j+1, k, mcCorrMin );

				binMigration_full->SetBinContent( i+1, j+1, k, binMigrationCorr );
				accCorrection_full->SetBinContent( i+1, j+1, k, accCorr );
				mcCorrection_full->SetBinContent( i+1, j+1, k, mcCorr );
				
				binMigration_p->SetBinError( i+1, j+1, k,  getUncertainty( matchBinPos, recBinPos ) );
				accCorrection_p->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, matchBinPos ) );
				mcCorrection_p->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, recBinPos ) );
				
				binMigration_m->SetBinError( i+1, j+1, k, getUncertainty( matchBinMin, recBinMin ) );
				accCorrection_m->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, matchBinMin) );
				mcCorrection_m->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, recBinMin) );

				binMigration_full->SetBinError( i+1, j+1, k, getUncertainty( binMigrationPos, binMigrationMin, getUncertainty( matchBinPos, recBinPos ), getUncertainty( matchBinMin, recBinMin  ) ) );
				accCorrection_full->SetBinError( i+1, j+1, k, getUncertainty( accCorrPos, accCorrMin, getUncertainty(genBinPos, matchBinPos ),getUncertainty(genBinMin, matchBinMin )  ) );
				mcCorrection_full->SetBinError( i+1, j+1, k, getUncertainty( mcCorrPos, mcCorrMin, getUncertainty(genBinPos, recBinPos ),getUncertainty(genBinMin, recBinMin )  ) );

			}
		}
	}
	binMigration_p->Write();
	accCorrection_p->Write();
	mcCorrection_p->Write();
				
	binMigration_m->Write();
	accCorrection_m->Write();
	mcCorrection_m->Write();

	binMigration_full->Write();
	accCorrection_full->Write();
	mcCorrection_full->Write();

	outFile->Close();
	//file_1->Close();
	//file_2->Close();
	
	return 1;
}

double getUncertainty( double num, double den ){

	return (num/den)*sqrt( pow( sqrt(num)/num , 2 ) + pow( sqrt(den)/den , 2 ) );

}

double getUncertainty( double num, double den, double num_unc, double den_unc ){

	return (num/den)*sqrt( pow( num_unc/num , 2 ) + pow( den_unc/den , 2 ) );

}
