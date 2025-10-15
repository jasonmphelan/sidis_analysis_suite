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
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	//TString in_name = argv[1];
    TString out_name = argv[1];
	TString in_name = argv[2];

	int nBinsXb = 14;
	int nBinsQ2 = 12;
	int nBinsZ = 14;

	double widthXb = (.66 - .1)/( (double) nBinsXb );
	double widthQ2 = (8. - 2.)/( (double) nBinsQ2 );
	double widthZ = (1. - .3)/( (double) nBinsZ );

        //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name, "RECREATE");
	
	TH1F* hM_rho[2][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho_bac[2][nBinsXb][nBinsQ2][nBinsZ];

	TFile * rho_norms = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms"+in_name+".root");
	TH3F * hNorms[2];
	hNorms[0] = (TH3F*)rho_norms->Get("hNorm_pip");
	hNorms[1] = (TH3F*)rho_norms->Get("hNorm_pim");

	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsZ; z++ ){
				TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*widthQ2, 2 + (q+1)*widthQ2, .3 + z*widthZ, .3 + (z + 1)*widthZ);
				
				hM_rho[0][x][q][z] = new TH1F( Form("hM_rho_pip_%i_%i_%i", x, q, z), "pip " + binLims, 50, 0, 2.5);
				hM_rho[1][x][q][z] = new TH1F( Form("hM_rho_pim_%i_%i_%i", x, q, z), "pim " + binLims, 50, 0, 2.5);
				hM_rho_bac[0][x][q][z] = new TH1F( Form("hM_rho_bac_pip_%i_%i_%i", x, q, z), "pip " + binLims, 50, 0, 2.5);
				hM_rho_bac[1][x][q][z] = new TH1F( Form("hM_rho_bac_pim_%i_%i_%i", x, q, z), "pim " + binLims, 50, 0, 2.5);
				
			}
		}
	}
	TChain * file_rec = new TChain("ePi");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2" + in_name + ".root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4" + in_name + ".root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6" + in_name + ".root");
	TTreeReader reader( file_rec);

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions();
	anal.loadMatchingFunctions3D();

	//Load input tree
    TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderArray<double> rhoWeight(reader_rec, "rhoWeight");

	//Fill histograms
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%10000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			if( !isGoodPion[i] ){continue;}
			//if( !isGoodPion[0] || !isGoodPion[1]){continue;}
			//if ( !anal.applyAcceptanceMatching(pi[0], 2) || !anal.applyAcceptanceMatching(pi[1], 2)) continue;

			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*nBinsXb);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*nBinsZ);
	

			//if( rhoWeight[i] > 10 ){ continue; }

			double scaling = hNorms[chargeIdx]->GetBinContent(this_bin_xB+1, this_bin_Q2+1, this_bin_Z+1);
			if( pi[0].getZ() + pi[1].getZ() < 0.5 || pi[0].getZ() + pi[1].getZ() > 1.1 ){continue;} 

			if( *Mx_2pi < 1.15 ){
				hM_rho[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) , rhoWeight[i] );
			}
			if( *Mx_2pi > 1.15 && *Mx_2pi < 1.35 ){
				hM_rho_bac[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho), rhoWeight[i]*scaling );
			}
			
		}
	}
	
	file_out->cd();
	
	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsZ; z++ ){
				hM_rho[0][x][q][z]->Write();
				hM_rho[1][x][q][z]->Write();
				hM_rho_bac[0][x][q][z]->Write();
				hM_rho_bac[1][x][q][z]->Write();
			}
		}
	}
	
	file_out->Close();

}
