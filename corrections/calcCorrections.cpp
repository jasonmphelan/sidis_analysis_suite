#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"


#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<clashit>+;
#pragma link C++ class vector<clashit>+;
#endif

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

const int nXbBins = 10;
const int nQ2Bins = 12;
const int nZBins = 14;

const double xB_min = .1;
const double xB_max = .6;
const double xB_width = (xB_max - xB_min)/( (double) nXbBins);

const double Q2_min = 2;
const double Q2_max = 8;
const double Q2_width = (Q2_max - Q2_min)/( (double) nQ2Bins);

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

        TTreeReaderValue<double> Q2_ptr(reader_rec, "Q2");
        TTreeReaderValue<double> xB_ptr(reader_rec, "xB");

        TTreeReaderArray<double> Z(reader_rec, "Z");
	TTreeReaderArray<int> charge(reader_rec, "charge");
        
	TTreeReaderValue<double> Q2_match_ptr(reader_rec, "Q2_gen");
        TTreeReaderValue<double> xB_match_ptr(reader_rec, "xB_gen");

        TTreeReaderArray<double> Z_match(reader_rec, "Z_gen");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");

	int event_total = recChain->GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}

                double Q2 = *Q2_ptr;
                double xB = *xB_ptr;
		double Q2_MC = *Q2_match_ptr;
		double xB_MC = *xB_match_ptr;
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
                int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);

                int this_bin_Q2_MC = (int)( ( (Q2_MC - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
                int this_bin_xB_MC = (int)( ( (xB_MC - xB_min)/(xB_max-xB_min) )*nXbBins);
		
	
		int pi_count = -1;
		for( auto Zval : Z ){
			pi_count++;
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion[pi_count]; }
			else if( matchType == 3 ){ matching = !isGoodPion3d[pi_count]; }
			else{ matching = false; }

			if( matching ){ continue; }

			int chargeIdx = (int)( charge[pi_count] < 1 );
		
			//Fill reco pions
			recHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( Zval );
			

			//Fill matched pions
			if( Z_match[pi_count] < .3 || Z_match[pi_count] > 1 ){ continue; }
			if( this_bin_Q2_MC < 0 || this_bin_Q2_MC >= nQ2Bins ){ continue; }
			if( this_bin_xB_MC < 0 || this_bin_xB_MC >= nXbBins ){ continue; }
			matchHists[this_bin_xB_MC][this_bin_Q2_MC][chargeIdx]->Fill( Z_match[pi_count] );
			
		}
	}
				


	
        TTreeReader reader_gen(genChain);

        TTreeReaderValue<double> Q2_MC_ptr(reader_gen, "Q2");
        TTreeReaderValue<double> xB_MC_ptr(reader_gen, "xB");

        TTreeReaderArray<double> Z_MC(reader_gen, "Z");
	TTreeReaderArray<int> charge_MC(reader_gen, "charge");
	
	TTreeReaderArray<bool> isGoodPion_MC(reader_gen, "isGoodPion");
	TTreeReaderArray<bool> isGoodPion3d_MC(reader_gen, "isGoodPion_3d");

	event_total = genChain->GetEntries();

	while (reader_gen.Next()) {
                int event_count = reader_gen.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<std::endl;
		}
		

                double Q2 = *Q2_MC_ptr;
                double xB = *xB_MC_ptr;
		
		int this_bin_Q2 = (int)( ( (Q2 - Q2_min)/(Q2_max-Q2_min) )*nQ2Bins);
                int this_bin_xB = (int)( ( (xB - xB_min)/(xB_max-xB_min) )*nXbBins);

		int pi_count = -1;
		for( auto Zval : Z_MC ){
			pi_count++;

	
			bool matching = true;

			if( matchType == 2 ){ matching = !isGoodPion_MC[pi_count]; }
			else if( matchType == 3 ){ matching = !isGoodPion3d_MC[pi_count]; }
			else{ matching = false; }

			if( matching ){ continue; }
			int chargeIdx = (int)( charge_MC[pi_count] < 1 );
			
			//Fill reco pions
			genHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( Zval );

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
				
				double binMigrationMin = matchBinMin/recBinMin;
				double accCorrMin = genBinMin/matchBinMin;				
	
				double accCorr = accCorrPos/accCorrMin;
				double binMigrationCorr = binMigrationPos/binMigrationMin;
				
				cout<<"acc : "<<accCorr<<std::endl;

				binMigration_p->SetBinContent( i+1, j+1, k,  binMigrationPos );
				accCorrection_p->SetBinContent( i+1, j+1, k,  accCorrPos );
				
				binMigration_m->SetBinContent( i+1, j+1, k, binMigrationMin );
				accCorrection_m->SetBinContent( i+1, j+1, k, accCorrMin );

				binMigration_full->SetBinContent( i+1, j+1, k, binMigrationCorr );
				accCorrection_full->SetBinContent( i+1, j+1, k, accCorr );
				
				binMigration_p->SetBinError( i+1, j+1, k,  getUncertainty( matchBinPos, recBinPos ) );
				accCorrection_p->SetBinError( i+1, j+1, k,  getUncertainty(genBinPos, matchBinPos ) );
				
				binMigration_m->SetBinError( i+1, j+1, k, getUncertainty( matchBinMin, recBinMin ) );
				accCorrection_m->SetBinError( i+1, j+1, k, getUncertainty( genBinMin, matchBinMin) );

				binMigration_full->SetBinError( i+1, j+1, k, getUncertainty( binMigrationPos, binMigrationMin, getUncertainty( matchBinPos, recBinPos ), getUncertainty( matchBinMin, recBinMin  ) ) );
				accCorrection_full->SetBinError( i+1, j+1, k, getUncertainty( accCorrPos, accCorrMin, getUncertainty(genBinPos, matchBinPos ),getUncertainty(genBinMin, matchBinMin )  ) );

			}
		}
	}
	binMigration_p->Write();
	accCorrection_p->Write();
				
	binMigration_m->Write();
	accCorrection_m->Write();

	binMigration_full->Write();
	accCorrection_full->Write();

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
