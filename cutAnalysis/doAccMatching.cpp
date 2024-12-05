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
#include "analyzer.h"
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
		cerr << "./code [out name] [in name]\n";
		return -1;
	}

	TString in_name = argv[2];
	TString outFileName = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

        TFile * file_rec = new TFile(in_name, "UPDATE");
	TFile * outFile = new TFile( outFileName, "RECREATE");

	TH1F * hTheta[6][7][2];	
	TH2F * hThetaP[6][2];	

	for( int sec = 0; sec < 6; sec++ ){
		hThetaP[sec][0] = new TH2F( Form("hTheta_P_sec_%i_pip", sec), "", 100, 1.25, 5, 100, 0, 40);
		hThetaP[sec][1] = new TH2F( Form("hTheta_P_sec_%i_pim", sec), "", 100, 1.25, 5, 100, 0, 40);
		for( int bin = 0; bin <= 6; bin++ ){
			hTheta[sec][bin][0] = new TH1F( Form("hTheta_sec_%i_bin_%i_pip", sec, bin), "", 1000, 0, 40);
			hTheta[sec][bin][1] = new TH1F( Form("hTheta_sec_%i_bin_%i_pim", sec, bin), "", 1000, 0, 40);
		}
	}

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

        TTreeReaderValue<electron> e(reader_rec, "e");

        TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion_no_acc");
	int event_total = reader_rec.GetEntries();

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%1000000 == 0){
			cout<<"Events Analyzed: "<<event_count<<" / "<<event_total<<endl;
		}
		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
                        int sector_i = pi[i].getDC_sector() - 1;
    			int chargeIdx = (int)(pi[i].getCharge() < 0);			
			
			if ( !isGoodPion[i] ) { continue; }
		
			double p = pi[i].get3Momentum().Mag();
			double theta = pi[i].get3Momentum().Theta()*rad_to_deg; 

			if(sector_i < 0){continue;}

			int this_bin_p = (int)( ( (p - p_min)/(p_max-p_min) )*NpBins);
			if(this_bin_p > 6){continue;}
			hTheta[sector_i][this_bin_p][chargeIdx]->Fill(theta);
			hThetaP[sector_i][chargeIdx]->Fill(p, theta);
		}        		
	}
		
	//Get Cut Value (cut out 10% percentile)
	
	double cutValsMin[2][6][7];
	double cutValsMax[2][6][6];
	
	outFile->cd();

	for( int sec = 0; sec < 6; sec++ ){
		hThetaP[sec][0]->Write();
		hThetaP[sec][1]->Write();
		for( int bin = 0; bin < 7; bin++ ){
			for( int idx = 0; idx < 2; idx++ ){
				cutValsMin[idx][sec][bin] = getThetaPct( .01, hTheta[sec][bin][idx] );
				if( bin > 0 ){
					cutValsMax[idx][sec][bin-1] = getThetaPct( .99, hTheta[sec][bin][idx] );
				}
			}
		}
	}

	//for( int sec = 0; sec < 6; sec++ ){
	//	cout<<"Cut Vals Sector "<<sec<<" :\n";
	//	for( int idx = 0; idx < 2; idx++ ){
	///		for( int bin = 0; bin < 7; bin++ ){
	//			cout<<" "<<cutValsMin[idx][sec][bin]<<"\t";
	///		}
	///		cout<<"\n"<<endl;
	//		
	//		for( int bin = 0; bin < 6; bin++ ){
	//			cout<<" "<<cutValsMax[idx][sec][bin]<<"\t";
	//		}
	//		cout<<"\n\n\n\n"<<endl;
	//	}
	//}
	

	
	outFile->cd();

	double binVals[7] = {1.5, 2.0, 2.5, 3, 3.5, 4, 4.5};
	double binVals2[6] = {2.0, 2.5, 3, 3.5, 4, 4.5};
	//double binVals2[5] = { 2.5, 3, 3.5, 4, 4.5};
	TGraph *  gThetaPMax[6][2];
	TGraph *  gThetaPMin[6][2];
	TF1 *fitMax[6][2];
	TF1 *fitMin[6][2];
	TString label[2] = {"pip", "pim"};

	for( int sec = 0; sec < 6; sec++ ){
		for( int idx = 0; idx < 2; idx++ ){
			fitMin[sec][idx] = new TF1(Form("min_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);
			fitMin[sec][idx]->SetParameters(6, 7); 
			fitMin[sec][idx]->SetParNames("a", "b"); 

			fitMax[sec][idx] = new TF1(Form("max_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);
			fitMin[sec][idx]->SetParameters(20, 10); 
			fitMin[sec][idx]->SetParNames("a", "b"); 

			gThetaPMin[sec][idx] = new TGraph( 7, binVals, cutValsMin[idx][sec] );
			gThetaPMax[sec][idx] = new TGraph( 6, binVals2, cutValsMax[idx][sec] );
			
			cout<<"Fit Max : "<<idx<<endl;
			gThetaPMax[sec][idx]->Fit(Form("max_%i_", sec)+label[idx]);
			cout<<"Fit Min : "<<idx<<endl;
			gThetaPMin[sec][idx]->Fit(Form("min_%i_", sec)+label[idx]);

			fitMin[sec][idx]->Write();
			fitMax[sec][idx]->Write();
		}
	}

	TCanvas * c1 = new TCanvas("c1", "c1");
	hThetaP[1][0]->Draw("COLZ");
	fitMax[1][0]->Draw("SAME");
	fitMin[1][0]->Draw("SAME");
	
	fitMax[1][1]->SetLineColor(kGreen);
	fitMin[1][1]->SetLineColor(kGreen);
	
	fitMax[1][1]->Draw("SAME");
	fitMin[1][1]->Draw("SAME");

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
