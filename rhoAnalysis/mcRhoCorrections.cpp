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

const int nXbBins = bins_xB/14;
const int nQ2Bins = bins_Q2/12;
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
	
	TH1F * rhoHists[nXbBins][nQ2Bins][2];
	TH1F * allHists[nXbBins][nQ2Bins][2];
	
	TH3F * rhoCorrection_p = new TH3F( "hRhoCorrectionP", "hRhoCorrectionP", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * rhoCorrection_m = new TH3F( "hRhoCorrectionM", "hRhoCorrectionM", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	TH3F * rhoCorrection_full = new TH3F( "hRhoCorrection", "hRhoCorrection", nXbBins, xB_min, xB_max, nQ2Bins, Q2_min, Q2_max, nZBins, .3, 1); 
	
	
	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			
			rhoHists[i][j][0] = new TH1F( (TString)"rhoHist_P_"+Form("_%i_%i", i, j), Form("recHist_P_%i_%i", i, j), 14, .3, 1); 
			allHists[i][j][0] = new TH1F((TString)"allHist_P_"+Form("_%i_%i", i, j), Form("matchHist_P_%i_%i", i, j), 14, .3, 1); 
			
			rhoHists[i][j][1] = new TH1F( (TString)"rhoHist_M_"+Form("_%i_%i", i, j), Form("recHist_M_%i_%i", i, j), 14, .3, 1); 
			allHists[i][j][1] = new TH1F( (TString)"allHist_M_"+Form("_%i_%i", i, j), Form("matchHist_M_%i_%i", i, j), 14, .3, 1); 			
		}
	}

	TH1F * hMx = new TH1F("mx", "", 100, 0, 3);
	TH1F * hPx = new TH1F("px", "", 100, 0, 5);
	
	analyzer anal(0, -1);
	anal.setAnalyzerLevel(1);
	anal.loadMatchingFunctions();
	anal.loadMatchingFunctions3D();	
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));
	
	TTreeReader reader_rec(recChain);
	
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<bool> isGoodGenPion(reader_rec, "isGoodGenPion");

	TTreeReaderArray<bool> isGoodPion3d(reader_rec, "isGoodPion_3d");
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<genElectron> e_MC(reader_rec, "e_gen");
	TTreeReaderArray<pion> pi_vec(reader_rec, "pi");
	TTreeReaderArray<genPion> pi_match(reader_rec, "pi_gen");
	TTreeReaderArray<TLorentzVector> parent_p( reader_rec, "parent_p");
	
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
			if( (matchType==2 && isGoodPion[pi_count]) ){
				
				double Mx_2pi = (p_rest + e_MC->getQ() - parent_p[pi_count]).Mag2();
				double Px_2pi = (e_MC->getQ().Vect() - parent_p[pi_count].Vect()).Mag();
				hMx->Fill(sqrt(Mx_2pi));
				hPx->Fill(Px_2pi);
				//Fill reco pions
				allHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );
				if( (Mx_2pi < 0.1 || Mx_2pi > .95) && abs(pi.getParent()) != 211 ){
					rhoHists[this_bin_xB][this_bin_Q2][chargeIdx]->Fill( pi.getZ() );	
				}
			}
			
		}
	}
	
	outFile->cd();
	
	hMx->Write();
	hPx->Write();
	for( int i = 0; i < nXbBins; i++ ){
		for( int j = 0; j < nQ2Bins; j++ ){
			rhoHists[i][j][0]->Write();
			rhoHists[i][j][1]->Write();
			allHists[i][j][0]->Write();
			allHists[i][j][1]->Write();
			
			for( int k = 1; k <= nZBins; k++ ){
				
				double recBinPos = rhoHists[i][j][0]->GetBinContent(k);
				double recBinMin = rhoHists[i][j][1]->GetBinContent(k);
				
				double allBinPos = allHists[i][j][0]->GetBinContent(k);
				double allBinMin = allHists[i][j][1]->GetBinContent(k);
				
				double rhoCorrectionPos = recBinPos/allBinPos;
				double rhoCorrectionMin = recBinMin/allBinMin;
			
				double rhoCorrectionCorr = rhoCorrectionPos/rhoCorrectionMin;
				
				
				
				
				
				rhoCorrection_p->SetBinContent( i+1, j+1, k,  rhoCorrectionPos );
			
				rhoCorrection_m->SetBinContent( i+1, j+1, k, rhoCorrectionMin );
				
				rhoCorrection_full->SetBinContent( i+1, j+1, k, rhoCorrectionCorr );
				
				rhoCorrection_p->SetBinError( i+1, j+1, k,  getUncertainty( allBinPos, recBinPos ) );
				
				rhoCorrection_m->SetBinError( i+1, j+1, k, getUncertainty( allBinMin, recBinMin ) );
				
				
				rhoCorrection_full->SetBinError( i+1, j+1, k, getUncertainty( rhoCorrectionPos, rhoCorrectionMin, getUncertainty( allBinPos, recBinPos ), getUncertainty( allBinMin, recBinMin  ) ) );
				
			}
		}
	}
	rhoCorrection_p->Write();

	
	rhoCorrection_m->Write();

	
	rhoCorrection_full->Write();

	
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
