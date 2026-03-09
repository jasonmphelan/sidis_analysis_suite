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
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "constants.h"
#include "cut_values.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

using namespace cutVals;
using namespace constants;


int main( int argc, char** argv){

	if( argc <2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] <<"\n";

    TString out_name = argv[1];
	
	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions();
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));

    //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name, "RECREATE");
	
	TH1F* hWeights_pi[4][2][5];
	TH1F* hWeights_e[4][2][10];
	TString charge[2] = {"pip", "pim"};
	TString type[4] = {"event", "sym", "event_scaled", "sym_scaled"};
	for( int t = 0; t < 4; t++ ){
		for( int ch = 0; ch < 2; ch++ ){
			for( int p = 0; p < 10; p++ ){
					int nBins = (t < 2) ? 4000 : 2000;
					if( p < 5 ) hWeights_pi[t][ch][p] = new TH1F( Form("hWeights_pi_"+type[t]+"_"+charge[ch]+"_%i", p), "", nBins, 0, 1000);

					hWeights_e[t][ch][p] = new TH1F( Form("hWeights_e_"+type[t]+"_"+charge[ch]+"_%i", p), "", nBins, 0, 1000);
				
			}
		}
	}

	TH1F * hZ_sum = new TH1F( "hZ_sum", "", 100, 0.5, 1);

	TChain * file_rec = new TChain("ePi");
	file_rec->Add("../trees/final_skims/rho_skims/rotated_10.2_sym.root");
	file_rec->Add("../trees/final_skims/rho_skims/rotated_10.4_sym.root");
	file_rec->Add("../trees/final_skims/rho_skims/rotated_10.6_sym.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root");
	TTreeReader reader( file_rec);


	//Load input tree
        TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion"); //should be no acceptance
	TTreeReaderArray<double> rhoWeight(reader_rec, "rhoWeight");
	TTreeReaderArray<double> rhoWeight_sym(reader_rec, "rhoWeight_sym");
	

	//Fill histograms
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
		if( *Mx_2pi > 1.125){ continue;}
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			
			if( !isGoodPion[i] ){continue;}
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			double p_e = e->get3Momentum().Mag();
			double p_pi = pi[i].get3Momentum().Mag();

			if( !anal.applyAcceptanceMatching(pi[i], 2) ) continue;
			
			//if( rhoWeight[i] == 0 || rhoWeight[i] == 1 || rhoWeight[i] > 10 ){ continue; }
			if( p_pi > 5 || p_e > 10) continue;
			
			hWeights_e[0][chargeIdx][(int)(p_e)]->Fill( rhoWeight[i] );//, rhoWeight[i] );
			hWeights_pi[0][chargeIdx][(int)(p_pi)]->Fill( rhoWeight[i] );//, rhoWeight[i] );
			
			hWeights_e[1][chargeIdx][(int)(p_e)]->Fill( rhoWeight_sym[i] );//, rhoWeight[i] );
			hWeights_pi[1][chargeIdx][(int)(p_pi)]->Fill( rhoWeight_sym[i] );//, rhoWeight[i] );

			hWeights_pi[2][chargeIdx][1]->Fill( rhoWeight[i], rhoWeight[i] );

			hWeights_pi[3][chargeIdx][1]->Fill( rhoWeight_sym[i] , rhoWeight_sym[i] );
		}
		hZ_sum->Fill( pi[0].getZ() + pi[1].getZ() );
	}

	cout<<"writing hists\n";

	file_out->cd();
	for( int t =0; t < 4; t++ ){
		for( int ch = 0; ch < 2; ch++ ){
		
			for( int p = 0; p < 10; p++ ){
			
				if( p < 5 && p > 0 ) {
					cout<<"writing pi\n";
					hWeights_pi[t][ch][p]->Write();
				}
			
				if( p > 2 ) hWeights_e[t][ch][p]->Write();
			}
		}
	}
	hZ_sum->Write();
	file_out->Close();

}
