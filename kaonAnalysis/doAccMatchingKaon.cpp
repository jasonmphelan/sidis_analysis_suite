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
#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
#include "reader.h"
#include "analyzer.h"
#include "reader.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TGraph.h"

using namespace clas12;
using namespace constants;

using namespace cutVals; 

const int NpBins = 7;

const double p_min = 1.25;
const double p_max = 4.75;

double getThetaPct( double pct, TH1F * h);
double fitF( double *p, double *par);

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [out name] [in type (pi2k = 0, k2pi = 1)]\n";
		return -1;
	}
	int cType = atoi(argv[2]);

	TChain * file_rec = new TChain("ePi");
	if( cType == 0 ){
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root");
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root");
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root");
	}
	if( cType == 1 ){
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.2/final_skim.root");
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.4/final_skim.root");
		file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/kaons_10.6/final_skim.root");
	}
	TTreeReader reader( file_rec);
	TString outFileName = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	TFile * outFile = new TFile( outFileName, "RECREATE");

	TH1F * hTheta[7][2];	
	TH2F * hThetaP[2];	

	hThetaP[0] = new TH2F( "hTheta_P_pip", "", 100, 1.25, 5, 100, 0, 40);
	hThetaP[1] = new TH2F( "hTheta_P_pim", "", 100, 1.25, 5, 100, 0, 40);
	for( int bin = 0; bin <= 6; bin++ ){
		hTheta[bin][0] = new TH1F( Form("hTheta_bin_%i_pip", bin), "", 1000, 0, 40);
		hTheta[bin][1] = new TH1F( Form("hTheta_bin_%i_pim", bin), "", 1000, 0, 40);
	}
	

	//Load input tree
        TTreeReader reader_rec(file_rec);

        TTreeReaderValue<electron> e(reader_rec, "e");

        TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");
	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<endl;
		}
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
    			int chargeIdx = (int)(pi[i].getCharge() < 0);			

			if ( pi[i].getBeta_rich() <.0001 ){ continue; }
			if ( !isGoodPion[i] ) { continue; }
			
			double p = pi[i].get3Momentum().Mag();
			double theta = pi[i].get3Momentum().Theta()*rad_to_deg; 


			int this_bin_p = (int)( ( (p - p_min)/(p_max-p_min) )*NpBins);
			if(this_bin_p > 6){continue;}
			hTheta[this_bin_p][chargeIdx]->Fill(theta);
			hThetaP[chargeIdx]->Fill(p, theta);
		}        		
	}
		
	//Get Cut Value (cut out 1% percentile)
	
	double cutValsMin_pim[7];
	double cutValsMax_pim[7]; //note, we skip first p bin for calc of max
	
	double cutValsMin_pip[5];
	double cutValsMax_pip[5]; //note, we skip first p bin for calc of max
	
	outFile->cd();

	hThetaP[0]->Write();
	hThetaP[1]->Write();
	for( int bin = 0; bin < 7; bin++ ){
		cutValsMin_pim[bin] = getThetaPct( .01, hTheta[bin][1] );
		
		cutValsMax_pim[bin] = getThetaPct( .99, hTheta[bin][1] );	
		
		if( bin >= 2 ){
			cutValsMin_pip[bin-2] = getThetaPct( .01, hTheta[bin][0] );
			cutValsMax_pip[bin-2] = getThetaPct( .99, hTheta[bin][0] );	
			cout<<cutValsMax_pip[bin-2]<<std::endl;
		}	
	}

	outFile->cd();

	double binVals_pim_max[7] = {1.5, 2.0, 2.5, 3, 3.5, 4, 4.5};
	double binVals_pim_min[7] = {1.5, 2.0, 2.5, 3, 3.5, 4, 4.5};
	double binVals_pip_max[5] = {2.5, 3, 3.5, 4, 4.5};
	double binVals_pip_min[5] = {2.5, 3, 3.5, 4, 4.5};
	//double binVals2[5] = { 2.5, 3, 3.5, 4, 4.5};
	TGraph *  gThetaPMax[2];
	TGraph *  gThetaPMin[2];
	TF1 *fitMax[6][2];
	TF1 *fitMin[6][2];
	TString label[2] = {"pip", "pim"};
	double lims[2] = {2.5, 1.25};	


	for( int sec = 0; sec < 6; sec++ ){
		for( int idx = 0; idx < 2; idx++ ){
			fitMin[sec][idx] = new TF1(Form("min_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);
			fitMin[sec][idx]->SetParameters(6, 0); 
			fitMin[sec][idx]->SetParNames("a", "b"); 

			fitMax[sec][idx] = new TF1(Form("max_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);
			fitMax[sec][idx]->SetParameters(20, 0); 
			fitMax[sec][idx]->SetParNames("a", "b"); 
		
			if(sec != 0){
				fitMin[sec][idx]->SetParameters(0, -999);
				fitMax[sec][idx]->SetParameters(999, 999);
				fitMin[sec][idx]->Write();
				fitMax[sec][idx]->Write();
			}
			else{
				if( idx == 1 ){
					gThetaPMin[idx] = new TGraph( 7, binVals_pim_min, cutValsMin_pim );
					gThetaPMax[idx] = new TGraph( 7, binVals_pim_max, cutValsMax_pim );
				}
				else{
					gThetaPMin[idx] = new TGraph( 5, binVals_pip_min , cutValsMin_pip );
					gThetaPMax[idx] = new TGraph( 5, binVals_pip_max, cutValsMax_pip );
				}

				gThetaPMax[idx]->Write( (TString)"gMax" + label[idx]);
				gThetaPMin[idx]->Write( (TString)"gMin" + label[idx]);
				
				gThetaPMax[idx]->Fit(Form("max_%i_", sec)+label[idx]);
				gThetaPMin[idx]->Fit(Form("min_%i_", sec)+label[idx]);

				fitMin[sec][idx]->Write();
				fitMax[sec][idx]->Write();
			}
		}
	}
	outFile->Close();
}

double getThetaPct( double pct, TH1F * h){
	double cutoff = h->Integral()*pct;
	
	double cumCount = 0;

	for ( int i = 1; i <= h->GetNbinsX(); i++ ){
		cumCount += h->GetBinContent(i);
		
		if( cumCount > cutoff ){
			return h->GetBinCenter(i);
		}
	}

	return h->GetBinCenter( h->GetNbinsX() );
}

double fitF( double *p, double *par){
	return par[0] + par[1]/p[0];
}
