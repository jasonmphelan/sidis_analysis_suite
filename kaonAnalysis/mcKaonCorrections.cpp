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
const int nZBins = 2*bins_Z;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );

int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Rec File] [Output File] [Matching]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2]  << "\n";

	TString inName_rec = argv[1];
	TString outName = argv[2];
	int matchType = atoi( argv[3] );

	TFile * outFile = new TFile( outName, "RECREATE");

	TFile * inFile_rec = new TFile( inName_rec );
	TTree * recChain = (TTree *)inFile_rec->Get("ePi");
	

	TH1F * recHists[4][nXbBins][nQ2Bins][2];
	TH1F * matchHists[4][nXbBins][nQ2Bins][2];

	TH3F * mcCorrection_p[4];
	TH3F * mcCorrection_m[4];
	TH3F * mcCorrection_full[4];

	for( int p = 0; p < 4; p++ ){

		TH3F * mcCorrection_p[p] = new TH3F( Form("hKaonCorrP_%i", p), "hMcCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		TH3F * mcCorrection_m[p] = new TH3F( Form("hKaonCorrM_%i",p), "hMcCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		TH3F * mcCorrection_full[p] = new TH3F( "hKaonCorr", "hMcCorrection", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	

		for( int i = 0; i < nXbBins; i++ ){
			for( int j = 0; j < nQ2Bins; j++ ){
			
				recHists[p][i][j][0] = new TH1F( (TString)"recHist_P_"+Form("_%i_%i", i, j), Form("recHist_P_%i_%i", i, j), 14, .3, 1); 
				matchHists[p][i][j][0] = new TH1F((TString)"matchHist_P_"+Form("_%i_%i", i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1); 
				
				recHists[p][i][j][1] = new TH1F( (TString)"recHist_M_"+Form("_%i_%i", i, j), Form("recHist_M_%i_%i", i, j), 14, .3, 1); 
				matchHists[p][i][j][1] = new TH1F( (TString)"matchHist_M_"+Form("_%i_%i", i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1); 
			
			}
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
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
                int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);

		int pi_count = -1;
		for( auto pi : pi_vec ){
			pi_count++;
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[pi_count]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[pi_count]; }
			else{ matching = false; }

			if( matching ){ continue; }

			int chargeIdx = (int)( pi.getCharge() < 1 );
			int this_bin_p = -1;
			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
					this_bin_p = j;
				}
			}
			if( this_bin_p<0 ){continue;}

			//Fill reco pions
			recHists[this_bin_p][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
			

			//Fill matched pions
			if( pi.getPID() == 211 ){
				matchHists[this_bin_p][this_bin_xB_MC][this_bin_Q2_MC][chargeIdx]->Fill( pi_match[pi_count].getZ() );
			}
		}
	}
				


	outFile->cd();


	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			recHists[i][j][0]->Write();
			recHists[i][j][1]->Write();
			matchHists[i][j][0]->Write();
			matchHists[i][j][1]->Write();
			
			for( int k = 1; k <= nZBins; k++ ){

				double recBinPos = recHists[i][j][0]->GetBinContent(k);
				double recBinMin = recHists[i][j][1]->GetBinContent(k);

				double matchBinPos = matchHists[i][j][0]->GetBinContent(k);
				double matchBinMin = matchHists[i][j][1]->GetBinContent(k);

				double binMigrationPos = matchBinPos/recBinPos;

				double binMigrationMin = matchBinMin/recBinMin;
				double accCorrMin = genBinMin/matchBinMin;				
				double mcCorrMin = genBinMin/recBinMin;
	
				double accCorr = accCorrPos/accCorrMin;
				double binMigrationCorr = binMigrationPos/binMigrationMin;
				double mcCorr = mcCorrPos/mcCorrMin;
				
								


				mcCorrection_p->SetBinContent( i+1, j+1, k,  mcCorrPos );
			
				mcCorrection_m->SetBinContent( i+1, j+1, k, mcCorrMin );

				mcCorrection_full->SetBinContent( i+1, j+1, k, mcCorr );
				
				mcCorrection_p->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, recBinPos ) );
				
				mcCorrection_m->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, recBinMin) );

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
	
	return 1;
}

double getUncertainty( double num, double den ){

	return (num/den)*sqrt( pow( sqrt(num)/num , 2 ) + pow( sqrt(den)/den , 2 ) );

}

double getUncertainty( double num, double den, double num_unc, double den_unc ){

	return (num/den)*sqrt( pow( num_unc/num , 2 ) + pow( den_unc/den , 2 ) );

}
