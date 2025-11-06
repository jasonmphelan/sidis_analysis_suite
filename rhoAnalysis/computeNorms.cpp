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
using std::isnan;
using namespace cutVals;
using namespace constants;


int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
    TString out_name = argv[2];
	double Mx_low = atof(argv[3]);
	double Mx_high = atof(argv[4]);
	double norm_bound = atof(argv[5]);

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));


	int nBinsXb = 14;
	int nBinsQ2 = 12;
	int nBinsZ = 14;

	double widthXb = (.66 - .1)/( (double) nBinsXb );
	double widthQ2 = (8. - 2.)/( (double) nBinsQ2 );
	double widthZ = (1. - .3)/( (double) nBinsZ );

        //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name + in_name + ".root", "RECREATE");
	TVector3 * bounds = new TVector3(Mx_low, Mx_high, norm_bound);

	TH1F* hMx_2pi[2][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho[2][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho_back[2][nBinsXb][nBinsQ2][nBinsZ];


	TH3F * hNorms[2];
	hNorms[0] = new TH3F( "hNorm_pip", "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);
	hNorms[1] = new TH3F( "hNorm_pim", "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);

	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsZ; z++ ){
				TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*widthQ2, 2 + (q+1)*widthQ2, .3 + z*widthZ, .3 + (z + 1)*widthZ);
				hMx_2pi[0][x][q][z] = new TH1F( Form("hMx_pip_%i_%i_%i", x, q, z), "pip " + binLims, 75, 0, 2.5);
				hMx_2pi[1][x][q][z] = new TH1F( Form("hMx_pim_%i_%i_%i", x, q, z), "pim " + binLims, 75, 0, 2.5);
				
				hM_rho[0][x][q][z] = new TH1F( Form("hM_rho_pip_%i_%i_%i", x, q, z), "pip " + binLims, 75, 0, 2.5);
				hM_rho[1][x][q][z] = new TH1F( Form("hM_rho_pim_%i_%i_%i", x, q, z), "pim " + binLims, 75, 0, 2.5);
				hM_rho_back[0][x][q][z] = new TH1F( Form("hM_rho_back_pip_%i_%i_%i", x, q, z), "pip " + binLims, 75, 0, 2.5);
				hM_rho_back[1][x][q][z] = new TH1F( Form("hM_rho_back_pim_%i_%i_%i", x, q, z), "pim " + binLims, 75, 0, 2.5);
			}
		}
	}
	TChain * file_rec = new TChain("ePi");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2" + in_name + ".root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4" + in_name + ".root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6" + in_name + ".root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root");


	//Load input tree
    TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion"); //should be no acceptance
	TTreeReaderArray<double> rhoWeight(reader_rec, "rhoWeight");

	//Fill histograms
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			if( !isGoodPion[i] ){continue;}
			//if( !isGoodPion[0] || !isGoodPion[1]){continue;}
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*nBinsXb);
			
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*nBinsZ);
	
			
			//bool matching = anal.applyAcceptanceMatching(pi[i], 2);
				//matching = isGoodPion[i]; }
			
			//if( !matching ){ continue; }


			//if( rhoWeight[i] <= 0 || rhoWeight[i] > 20 ){ continue; }
		
			hMx_2pi[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*Mx_2pi) , rhoWeight[i] );
			if( *Mx_2pi > bounds->X() && *Mx_2pi < bounds->Y() ) hM_rho_back[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) );//, rhoWeight[i] );
			if( *Mx_2pi < bounds->X() ) hM_rho[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) );//, rhoWeight[i] );

		}
	}

	
	TFile * out_fits = new TFile(out_name + "_M_rho.root", "RECREATE");
	out_fits->cd();

	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsZ; z++ ){				
			
				double bin_max= hM_rho[0][x][q][z]->FindBin( bounds->Z() );
				double back_num_pip = 0;
				double back_num_pim = 0;
				double rho_num_pip = 0;
				double rho_num_pim = 0;

				for( int bin = 1; bin <= bin_max; bin++ ){
					back_num_pip+=hM_rho_back[0][x][q][z]->GetBinContent(bin);
					back_num_pim+=hM_rho_back[1][x][q][z]->GetBinContent(bin);
					rho_num_pip+=hM_rho[0][x][q][z]->GetBinContent(bin);
					rho_num_pim+=hM_rho[1][x][q][z]->GetBinContent(bin);
				}


				double norm_pip = rho_num_pip/back_num_pip;//fit_dis_pip->Integral( .75, 1.15 )/fit_dis_pip->Integral(1.15, 1.45);
				double norm_pim = rho_num_pim/back_num_pim;//fit_dis_pim->Integral( .75, 1.15 )/fit_dis_pim->Integral(1.15, 1.45);
				if( norm_pip < 0 || isnan(norm_pip) || !isfinite(norm_pip) ) norm_pip = 0;
				if( norm_pim < 0 || isnan(norm_pim) || !isfinite(norm_pim) ) norm_pim = 0;
				
				if( norm_pip > 0 && !isnan(norm_pip) && isfinite(norm_pip) )cout<<"Norm pip "<<norm_pip<<std::endl;

				hNorms[0]->SetBinContent( x+1, q+1, z+1, norm_pip );
				hNorms[1]->SetBinContent( x+1, q+1, z+1, norm_pim );
				
				hMx_2pi[0][x][q][z]->Write();
				hMx_2pi[1][x][q][z]->Write();
			}
		}
	}
	
	out_fits->Close();
	
	file_out->cd();
	bounds->Write("bounds");
	hNorms[0]->Write();
	hNorms[1]->Write();
	file_out->Close();

}
