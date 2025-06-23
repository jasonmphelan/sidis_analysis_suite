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

const int nXbBins = bins_xB/2;
const int nQ2Bins = bins_Q2;
const int nZBins = 2*bins_Z;

double getUncertainty( double num, double den );
double getUncertainty( double num, double den, double num_unc, double den_unc );

int main( int argc, char** argv){

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Rec File] [Output File] [Matching] [pi2k (0), k2pi(1)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2]  << "\n";

	TString inName_rec = argv[1];
	TString outName = argv[2];
	int matchType = atoi( argv[3] );
	int corrType = atoi( argv[4] );
	int pid = 211;
	if( corrType == 1 ){
		pid = 321;
	}

	TFile * outFile = new TFile( outName, "RECREATE");

	TFile * inFile_rec = new TFile( inName_rec );
	TTree * recChain = (TTree *)inFile_rec->Get("ePi");
	

	TH1F * allHists[4][nXbBins][nQ2Bins][2];
	TH1F * piHists[4][nXbBins][nQ2Bins][2];

	TH3F * mcCorrection_p[4];
	TH3F * mcCorrection_m[4];
	TH3F * mcCorrection_full[4];

	for( int p = 0; p < 4; p++ ){

		mcCorrection_p[p] = new TH3F( Form("hKaonCorrP_%i", p), "hKaonCorrP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		mcCorrection_m[p] = new TH3F( Form("hKaonCorrM_%i",p), "hKaonCorrM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
		mcCorrection_full[p] = new TH3F( Form("hKaonCorr_%i", p), "hKaonCorr", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	

		for( int i = 0; i < nXbBins; i++ ){
			for( int j = 0; j < nQ2Bins; j++ ){
			
				allHists[p][i][j][0] = new TH1F( (TString)"recHist_P_"+Form("_%i_%i_%i", p, i, j), Form("recHist_P_%i_%i", i, j), nZBins, .3, 1); 
				piHists[p][i][j][0] = new TH1F((TString)"matchHist_P_"+Form("_%i_%i_%i", p, i, j), Form("matchHist_P_%i_%i", i, j), nZBins, .3, 1); 
				
				allHists[p][i][j][1] = new TH1F( (TString)"recHist_M_"+Form("_%i_%i_%i", p, i, j), Form("recHist_M_%i_%i", i, j), nZBins, .3, 1); 
				piHists[p][i][j][1] = new TH1F( (TString)"matchHist_M_"+Form("_%i_%i_%i", p, i, j), Form("matchHist_M_%i_%i", i, j), nZBins, .3, 1); 
			
			}
		}
	}

        TTreeReader reader_rec(recChain);

        
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderValue<electron> e(reader_rec, "e");

    TTreeReaderArray<pion> pi_vec(reader_rec, "pi");

	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");

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
			if( abs( pi.getPID_eb() ) != pid ){ continue; }
	
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[pi_count]; }
			else if( matchType == 3 ){ matching = !isGoodPion3d[pi_count]; }
			else{ matching = false; }

			if( matching ){ continue; }
			int chargeIdx = (int)( pi.getCharge() < 1 );
			int this_bin_p = -1;
			double p_pi = pi.get3Momentum().Mag();
			for( int j= 0; j < bins_p; j++ ){
				if( p_pi > p_bin_edges[j] && p_pi < p_bin_edges[j+1] ){
					this_bin_p = j;
				}
			}
			if( this_bin_p<0 ){continue;}

			//Fill reco pions

			allHists[this_bin_p][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
			

			//Fill matched pions
			if( abs(pi.getPID()) == 211 ){
				piHists[this_bin_p][this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
			}
		}
	}
				


	outFile->cd();

	for( int p = 0; p < 4; p++ ){
	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			allHists[p][i][j][0]->Write();
			allHists[p][i][j][1]->Write();
			piHists[p][i][j][0]->Write();
			piHists[p][i][j][1]->Write();
			
			for( int k = 1; k <= nZBins; k++ ){

				double allBinPos = allHists[p][i][j][0]->GetBinContent(k);
				double allBinMin = allHists[p][i][j][1]->GetBinContent(k);

				double piBinPos = piHists[p][i][j][0]->GetBinContent(k);
				double piBinMin = piHists[p][i][j][1]->GetBinContent(k);

				double mcCorrPos = piBinPos/allBinPos;	
				double mcCorrMin = piBinMin/allBinMin;
		
				double mcCorr = mcCorrPos/mcCorrMin;

				mcCorrection_p[p]->SetBinContent( i+1, j+1, k,  mcCorrPos );
			
				mcCorrection_m[p]->SetBinContent( i+1, j+1, k, mcCorrMin );

				mcCorrection_full[p]->SetBinContent( i+1, j+1, k, mcCorr );
				
				mcCorrection_p[p]->SetBinError( i+1, j+1, k,  getUncertainty(piBinPos, allBinPos ) );
				
				mcCorrection_m[p]->SetBinError( i+1, j+1, k, getUncertainty( piBinMin, allBinMin) );

				mcCorrection_full[p]->SetBinError( i+1, j+1, k, getUncertainty( mcCorrPos, mcCorrMin, getUncertainty(piBinPos, allBinPos ) ,getUncertainty(piBinMin, allBinMin )   ) );

			}
		}
	}

	mcCorrection_p[p]->Write();
	mcCorrection_m[p]->Write();
	mcCorrection_full[p]->Write();
}
	outFile->Close();
	
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

