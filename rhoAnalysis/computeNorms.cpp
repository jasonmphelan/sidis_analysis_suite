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

	TString in_name = argv[1];
       	TString out_name = argv[2];

	int nBinsXb = 14;
	int nBinsQ2 = 12;
	int nBinsZ = 14;

	double widthXb = (.66 - .1)/( (double) nBinsXb );
	double widthQ2 = (8. - 2.)/( (double) nBinsQ2 );
	double widthZ = (1. - .3)/( (double) nBinsZ );

        //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name, "RECREATE");
	
	TH1F* hMx_2pi[2][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho[2][nBinsXb][nBinsQ2][nBinsZ];

	TH3F * hNorms[2];
	hNorms[0] = new TH3F( "hNorm_pip", "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);
	hNorms[1] = new TH3F( "hNorm_pim", "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);

	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsQ2; z++ ){
				TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*widthQ2, 2 + (q+1)*widthQ2, .3 + z*widthZ, .3 + (z + 1)*widthZ);
				hMx_2pi[0][x][q][z] = new TH1F( Form("hMx_pip_%i_%i_%i", x, q, z), "pip " + binLims, 100, 0, 2.5);
				hMx_2pi[1][x][q][z] = new TH1F( Form("hMx_pim_%i_%i_%i", x, q, z), "pim " + binLims, 100, 0, 2.5);
				
				hM_rho[0][x][q][z] = new TH1F( Form("hM_rho_pip_%i_%i_%i", x, q, z), "pip " + binLims, 100, 0, 2.5);
				hM_rho[1][x][q][z] = new TH1F( Form("hM_rho_pim_%i_%i_%i", x, q, z), "pim " + binLims, 100, 0, 2.5);
			}
		}
	}
	TChain * file_rec = new TChain("ePi");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2.root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4.root");
	file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6.root");
	TTreeReader reader( file_rec);


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
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*nBinsXb);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*nBinsZ);
	

			if( rhoWeight[i] <= 1 || rhoWeight[i] > 10 ){ continue; }

			hMx_2pi[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( *Mx_2pi, rhoWeight[i] );
			hM_rho[chargeIdx][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( *M_rho, rhoWeight[i] );
			
		}
	}

	//Fit histograms
	TString fitString = "[0]*exp( -([1] - x)*([1] - x)/(2.*[2]*[2] ) ) + [3]*exp( -([4] - x)*([4] - x)/(2.*[5]*[5] ) )";
	TF1 * fitFunc[2][nBinsXb][nBinsQ2][nBinsZ];
	
	TCanvas canvas("canvas");
	canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_fits.pdf[");
	canvas.Clear();
	
	for( int x = 0; x < nBinsXb; x++ ){
		for( int q = 0; q < nBinsQ2; q++ ){
			for( int z = 0; z <  nBinsQ2; z++ ){
				//TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*width, 2 + (q+1)*width, .3 + z*widthZ, .3 + (z + 1)widthZ);
				fitFunc[0][x][q][z] = new TF1( Form("fMx_pip_%i_%i_%i", x, q, z), fitString, 0, 2.5);
				fitFunc[1][x][q][z] = new TF1( Form("fMx_pim_%i_%i_%i", x, q, z), fitString, 0, 2.5);
						
				fitFunc[0][x][q][z]->SetParameters( hMx_2pi[0][x][q][z]->GetMaximum(),
					       				hMx_2pi[0][x][q][z]->GetMean(),
								       	hMx_2pi[0][x][q][z]->GetRMS(), 
									hMx_2pi[0][x][q][z]->GetMaximum()/2., .9, .05); 	
				fitFunc[1][x][q][z]->SetParameters( hMx_2pi[1][x][q][z]->GetMaximum(),
					       				hMx_2pi[1][x][q][z]->GetMean(), 
									hMx_2pi[1][x][q][z]->GetRMS(), 
									hMx_2pi[1][x][q][z]->GetMaximum()/2., .9, .05); 	
				//Fit and print missing mass plot

				hMx_2pi[0][x][q][z]->Fit( fitFunc[0][x][q][z] );
				if( hMx_2pi[0][x][q][z]->Integral() > 0 ){
					hMx_2pi[0][x][q][z]->Draw();
					canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_fits.pdf[");
					canvas.Clear();
				}
				hMx_2pi[1][x][q][z]->Fit( fitFunc[1][x][q][z] );
				if( hMx_2pi[1][x][q][z]->Integral() > 0 ){
					hMx_2pi[1][x][q][z]->Draw();
					canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_fits.pdf[");
					canvas.Clear();
				}

				canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_fits.pdf[");
				canvas.Clear();
				
				//extract norm
				TF1 * fit_temp_pip = new TF1("temp_pip", "gaus", .75, 1.15);
				TF1 * fit_temp_pim = new TF1("temp_pim", "gaus", .75, 1.45);
	
				fit_temp_pip->SetParameters( fitFunc[0][x][q][z]->GetParameter(0), fitFunc[0][x][q][z]->GetParameter(1), fitFunc[0][x][q][z]->GetParameter(2) );
				fit_temp_pim->SetParameters( fitFunc[1][x][q][z]->GetParameter(0), fitFunc[1][x][q][z]->GetParameter(1), fitFunc[1][x][q][z]->GetParameter(2) );
			
				double norm_pip = fit_temp_pip->Integral( .75, 1.15 )/fit_temp_pip->Integral(1.15, 1.45);
				double norm_pim = fit_temp_pim->Integral( .75, 1.15 )/fit_temp_pim->Integral(1.15, 1.45);
				if( norm_pip < 0 || isnan(norm_pip) || !isfinite(norm_pip) ) norm_pip = 0;
				if( norm_pim < 0 || isnan(norm_pim) || !isfinite(norm_pim) ) norm_pim = 0;
				hNorms[0]->SetBinContent( x+1, q+1, z+1, norm_pip );
				hNorms[1]->SetBinContent( x+1, q+1, z+1, norm_pim );
				
				//print results
				cout<<"Bin  : "<<x<< " , "<<q<<" , "<<z<<endl;
				cout<<"Pip norm : "<<hNorms[0]->GetBinContent( x+1, q+1, z+1 )<<endl;
				cout<<"Pim norm : "<<norm_pip<<endl;
			}
		}
	}
	canvas.Print("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/rho_fits.pdf]");
	
	
	file_out->cd();
	hNorms[0]->Write();
	hNorms[1]->Write();
	file_out->Close();

}
