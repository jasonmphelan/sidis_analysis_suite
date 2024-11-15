#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
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
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"

#define HIST_PATH _HIST
#define DATA_PATH _CORR

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File (no extension)]\n";
		cerr << "[Sector cut (optional, 0 if all)] [Theta Cut (optional)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

       	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
	
	cout<<"Creating Histograms\n";


	TH3F * hPhaseSpace[2];
	TH2F * hCorrVal[2][5];

	TString data_type[2] = {"pip", "pim"};
	TString pi_type[2] = {"#pi +", "#pi -"};

	for( int i = 0; i < 2; i++ ){
		hPhaseSpace[i] = new TH3F( Form("hPhaseSpace_%s", data_type[i].Data()), 
						";xB;Q2;z", 5, xB_min, xB_max, 
						50, Q2_min, Q2_max, 
						50, Z_min, Z_max);
		for( int j = 0; j < 5; j++ ){
			hCorrVal[i][j] = new TH2F( Form("hCorrVal_%s_%i", data_type[i].Data(), j ), 
							";z;Q^{2}", 500, Z_min, Z_max,
							500, Q2_min, Q2_max );
			hCorrVal[i][j]->SetTitle( Form("%s, x_{B} = %f", pi_type[i].Data(), ( (xB_min + .1*j) + (xB_min + .1*(j+1)) )/2. ) );
		}
	}


	cout<<"Beginning Event Loop\n";

	TFile * file = new TFile(in_name);
	//TChain * file = new TChain("ePi");
	//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.2/final_skim.root");
//	file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root");
//	file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root");
	TTreeReader reader("ePi", file);

	TTreeReaderValue<electron> e(reader, "e");

	TTreeReaderArray<pion> pi(reader, "pi");
	
	TTreeReaderArray<bool> isGoodPion_vec(reader, "isGoodPion");


	int event_count = 0;
	while (reader.Next()) {
                if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
		if( event_count == 10000000 ) break;       
		event_count++;
                double Q2 = e->getQ2();
                double xB = e->getXb();


		int chargeIdx = 0;
		
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			if( !isGoodPion_vec[i]) {continue;}
			chargeIdx = (int)(pi[i].getCharge() < 1);
			double Z = pi[i].getZ();

			hPhaseSpace[chargeIdx]->Fill(xB, Q2, Z);
		}
	}

	//For x value, get Corr values for filled bins
	correctionTools corr(0);
	cout<<"Load corrector\n";
	corr.loadParameters();

	cout<<"BeginPlotting\n";

	for( int x_bin = 0; x_bin < 5; x_bin++ ){
		double x_val = 	( (xB_min + .1*x_bin) + (xB_min + .1*(x_bin+1)) )/2.;
		
		for( int q_bin = 0; q_bin < 500; q_bin++ ){
			for( int z_bin = 0; z_bin < 500; z_bin++ ){
				for( int charge = 0; charge < 2; charge++ ){
					//Get kinematic values
					double q2_val = hCorrVal[charge][x_bin]->GetYaxis()->GetBinCenter( q_bin+1 );
					double z_val = hCorrVal[charge][x_bin]->GetXaxis()->GetBinCenter( z_bin+1 );
				
					//Check if bin is filled by data
					double q2_dat_bin = hPhaseSpace[charge]->GetYaxis()->FindBin( q2_val );
					double z_dat_bin = hPhaseSpace[charge]->GetZaxis()->FindBin( z_val );

					double bin_val = hPhaseSpace[charge]->GetBinContent( x_bin + 1, q2_dat_bin, z_dat_bin );
					if( bin_val < 1 ){ continue; }
					
					corr.setKinematics( x_val, q2_val, z_val, 0 );
					double correction = corr.getCorrectionFactor( 1, charge);
					cout<<"Correction value : "<<correction<<std::endl;
					hCorrVal[charge][x_bin]->Fill( z_val, q2_val, correction );
				}
			}
		}
	}




	//file->Close();

	outFile->cd();
	
	cout<<"BEGIN WRITING HISTOGRAMS\n";

	TCanvas canvas("canvas");
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf[");
	canvas.Clear();

	for( int j = 0; j < 5; j++ ){
		for( int i = 0; i < 2; i++ ){//Bin by charge
			hCorrVal[i][j]->Write();
			hCorrVal[i][j]->SetTitleSize(10);
			hCorrVal[i][j]->Draw("COLZ");
			canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf");
			canvas.Clear();
		}
			
	
	}
	canvas.Print((TString) HIST_PATH + "/" + out_name + ".pdf]");
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}
