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
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;


int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];

        TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name, "RECREATE");
	TH1F * hQ2[2];
	TH1F * hXb[2];
	TH1F * hZ[2];
	TH1F * hPt[2];
	TH1F * hMx_rho[2];
	TH1F * hMx_2pi[2];
	
	hQ2[0] = new TH1F("hQ2_all", ";Q^{2};", 50, Q2_min, Q2_max);
	hQ2[1] = new TH1F("hQ2_sideband", ";Q^{2};", 50, Q2_min, Q2_max);
	
	hMx_2pi[0] = new TH1F("hMx_2pi_all", ";Q^{2};", 50, .2, .9);
	hMx_2pi[1] = new TH1F("hMx_2pi_sideband", ";Q^{2};", 50, .2, .9);
	
	hMx_rho[0] = new TH1F("hMx_rho_all", ";Q^{2};", 50, .6, 1.4);
	hMx_rho[1] = new TH1F("hMx_rho_sideband", ";Q^{2};", 50, .6, 1.4);
	
	hXb[0] = new TH1F("hXb_all", ";x_{B};", 50, xB_min, xB_max);
	hXb[1] = new TH1F("hXb_sideband", ";x_{B};", 50, xB_min, xB_max);

	hZ[0] = new TH1F("hZ_all", ";z;", 50, Z_min, Z_max);
	hZ[1] = new TH1F("hZ_sideband", ";z;", 50, Z_min, Z_max);
	
	hPt[0] = new TH1F("hPt_all", ";p_{T};", 50, 0, 2);
	hPt[1] = new TH1F("hPt_sideband", ";p_{T};", 50, 0, 2);

	//Load input tree
        TTreeReader reader_rec("ePi", file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");

	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%500000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		hMx_2pi[0]->Fill(*Mx_2pi);
		hMx_rho[0]->Fill(*M_rho);
		
		if((*Mx_2pi > 1.1 && *Mx_2pi < 1.3) &&
		((*M_rho > .4 && *M_rho < .6)||
		(*M_rho > .85 && *M_rho < 1)) ){
			hMx_2pi[1]->Fill(*Mx_2pi);
			hMx_rho[1]->Fill(*M_rho);
		}
	
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			if( !isGoodPion[i] ){continue;}
			if( *Mx_2pi < 1.3 && *M_rho > .4 &&  *M_rho < 1){
				double z = (.77*.77 + (pi[0].get3Momentum() + pi[1].get3Momentum()).Mag2())/e->getOmega();			

				hQ2[0]->Fill( e->getQ2() );
				hXb[0]->Fill( e->getXb() );
				hZ[0]->Fill( z );
				hPt[0]->Fill( pi[i].getPi_q().Pt() );

			
			}

			if((*Mx_2pi > 1.1 && *Mx_2pi < 1.3) &&
			((*M_rho > .4 && *M_rho < .6)||
			(*M_rho > .85 && *M_rho < 1)) ){
				double z = (.77*.77 + (pi[0].get3Momentum() + pi[1].get3Momentum()).Mag2())/e->getOmega();			
			
				hQ2[1]->Fill( e->getQ2() );
				hXb[1]->Fill( e->getXb() );
				hZ[1]->Fill( z );
				hPt[1]->Fill( pi[i].getPi_q().Pt() );
			}
		}
	}

	file_out->cd();
	hQ2[0]->Write();
	hXb[0]->Write();
	hZ[0]->Write();
	hPt[0]->Write();
	
	hQ2[1]->Write();
	hXb[1]->Write();
	hZ[1]->Write();
	hPt[1]->Write();

	hMx_2pi[0]->Write();
	hMx_2pi[1]->Write();
	hMx_rho[0]->Write();
	hMx_rho[1]->Write();
	file_out->Close();

}
