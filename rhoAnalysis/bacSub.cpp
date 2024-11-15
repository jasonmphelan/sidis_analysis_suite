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

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	//TString in_name = argv[1];
       	TString out_name = argv[2];

	int nBinsXb = 14;
	int nBinsQ2 = 12;
	int nBinsZ = 14;

	double widthXb = (.66 - .1)/( (double) nBinsXb );
	double widthQ2 = (8. - 2.)/( (double) nBinsQ2 );
	double widthZ = (1. - .3)/( (double) nBinsZ );

        //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name, "RECREATE");
	
	TH1F* hM_rho_bac[2][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho[2][nBinsXb][nBinsQ2][nBinsZ];

	TH3F * hNorms[2];
	hNorms[0] = new TH3F( "hNorm_pip", "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);
	hNorms[1] = new TH3F( "hNorm_pim", "", nBinsXb, .1, .6, nBinsQ2, 2, 8, nBinsZ, .3, 1);

	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsQ2; z++ ){
				TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*widthQ2, 2 + (q+1)*widthQ2, .3 + z*widthZ, .3 + (z + 1)*widthZ);
				hM_rho[0][x][q][z] = new TH1F( Form("hM_rho_pip_%i_%i_%i", x, q, z), "pip " + binLims, 75, 0, 2.5);
				hM_rho[1][x][q][z] = new TH1F( Form("hM_rho_pim_%i_%i_%i", x, q, z), "pim " + binLims, 75, 0, 2.5);
				
				hM_rho_bac[0][x][q][z] = new TH1F( Form("hM_rho_bac_pip_%i_%i_%i", x, q, z), "pip " + binLims, 75, 0, 2.5);
				hM_rho_bac[1][x][q][z] = new TH1F( Form("hM_rho_bac_pim_%i_%i_%i", x, q, z), "pim " + binLims, 75, 0, 2.5);
			}
		}
	}


	//Load input tree
	TChain * file = new TChain("ePi");
	file->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_rotations/data_rotated_10.2.root");
	file->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_rotations/data_rotated_10.4.root");
	file->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_rotations/data_rotated_10.6.root");
	//TTreeReader reader( file);
        TTreeReader reader_rec(file);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	TTreeReaderValue<double> rhoWeight(reader_rec, "rhoWeight");

	//Fill histograms
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%10000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			if( !isGoodPion[i] ){continue;}
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*nBinsXb);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*nBinsZ);
		
			if( *rhoWeight <= 1 || *rhoWeight > 10 ){ continue; }

			if( *Mx_2pi > 1.15 && *Mx_2pi < 1.45 ){
				hM_rho_bac[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( *M_rho, *rhoWeight );
			}
			if( *Mx_2pi < 1.15 ){
				hM_rho[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( *M_rho, *rhoWeight );
			}
		}
	}

	TCanvas canvas("canvas");
	canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_bac.pdf[");
	canvas.Clear();

	TFile * normFile = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/build/test.root");
	TH3F * norms = (TH3F *)normFile->Get("hNorm_pip");
	cout<<"Do background subtraction"<<endl;
	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsQ2; z++ ){
				//TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*width, 2 + (q+1)*width, .3 + z*widthZ, .3 + (z + 1)widthZ);
				
				hM_rho_bac[0][x][q][z]->Scale( norms->GetBinContent( x+1, q+1, z+1 ) );

				hM_rho[0][x][q][z]->Draw("hist");
				hM_rho[0][x][q][z]->SetLineColor(kBlack);
				hM_rho_bac[0][x][q][z]->SetLineColor(kRed);
				hM_rho[0][x][q][z]->Add( hM_rho_bac[0][x][q][z], -1 );
				hM_rho[0][x][q][z]->SetLineColor(kAzure);
				hM_rho_bac[0][x][q][z]->Draw("same hist");
				hM_rho[0][x][q][z]->Draw("same hist");
				
				canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_bac.pdf");
				canvas.Clear();
				

			}
		}
	}
	canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_bac.pdf]");
	
	
	//file_out->cd();
	//file_out->Close();

}
