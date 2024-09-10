
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
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
#include "cut_values.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "reader.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;
using namespace cutVals; 

double getThetaPct( double pct, TH1D * h);
double fitF( double *p, double *par);

int main( int argc, char** argv){

	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input Path] [Output File] [# of input files] [File Type] [Beam Energy]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] << "\n";

	TString in_name = argv[1];
       	TString out_name = argv[2];
        TFile * inFile = new TFile(in_name);

	TFile * outFile = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/data/acceptanceMap.root", "RECREATE");

	TH2F * hThetaPhi[3][6][10]; //2D distribution
	TGraph * gThetaPhi[3][6][10]; //Graph of 99percentile
	TF1 * fThetaPhi[3][6][10]; //Fit functions
	TVector3 * fitBounds[3][6][10];

	TString momentum_functions = "[0]*(x - [1])*( x-[1] ) + [2]";
	//TString momentum_functions = "[0]*( -(x - [1])*( x-[1] ) )/( (x - [1])*( x-[1] ) - 1. ) + [2]";
	TString particle[3] = {"e","pip", "pim"};

	//loop through sectors, momentum, and particles	
	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < 10; bin++ ){
			for( int p = 0; p < 3; p++ ){
				//using fewer bins for pions
				if(p > 0 && bin >4){ continue; }

				hThetaPhi[p][sec][bin] = (TH2F *)inFile->Get( Form("hThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[p].Data()));
				//Set initial guesses
				double histMean = hThetaPhi[p][sec][bin]->GetMean(1);
				double histStd = hThetaPhi[p][sec][bin]->GetRMS(1);
				double histMeanY = hThetaPhi[p][sec][bin]->GetMean(2);
				double histStdY = hThetaPhi[p][sec][bin]->GetRMS(2);
				double histMin = histMeanY - 3*histStdY;
				

				fThetaPhi[p][sec][bin] = new TF1(Form("fThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[p].Data() ), 
								momentum_functions, histMean - 3*histStd, histMean + (3*histStd));
				fThetaPhi[p][sec][bin]->SetParameters( 2*histStd, histMean, histMin );
				fitBounds[p][sec][bin] = new TVector3( 0, 0, 0 );
			
			}

		}
			
	}
		
	//Get Cut Value (cut out 1% percentile)
		

	outFile->cd();


	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < 10; bin++ ){
			for( int i = 0; i < 3; i++ ){
				if( i > 0 && bin > 4 ){continue;}
				cout<<"NEW BIN-----------------------\n";
				int nPhiBins = hThetaPhi[i][sec][bin]->GetXaxis()->GetNbins();
				double binCenter[1000];
				double binVal[1000];
				int nPoints = 0;

				bool firstGoodBin = false;
				bool lastGoodBin = false;

				for( int phiBin = 0; phiBin < nPhiBins; phiBin++ ){
					TH1D * temp = new TH1D( "temp", "", hThetaPhi[i][sec][bin]->GetYaxis()->GetNbins(), 0, 40 );
					temp = (TH1D *)hThetaPhi[i][sec][bin]->ProjectionY( "temp", phiBin ,  phiBin   );
					if( temp->Integral() < 500 ){
						//set upper edge
						if( firstGoodBin && !lastGoodBin ){
							fitBounds[i][sec][bin]->SetY( hThetaPhi[i][sec][bin]->GetXaxis()->GetBinUpEdge( phiBin ) );
							lastGoodBin = true;
						}

						delete temp;
						continue;
					}
					//set lower edge
					if( !firstGoodBin ){
						fitBounds[i][sec][bin]->SetX( hThetaPhi[i][sec][bin]->GetXaxis()->GetBinLowEdge( phiBin ) );
						firstGoodBin = true;
					}
					
					//get data point based on percentile
					binCenter[nPoints] =  hThetaPhi[i][sec][bin]->GetXaxis()->GetBinCenter( phiBin )  ; 
					binVal[nPoints] = getThetaPct( .01, temp ) ;
					nPoints++;
			
					//set maximum
					double temp_max =  getThetaPct( 1, temp ) ;
					if( temp_max >  fitBounds[i][sec][bin]->Z() ){
						fitBounds[i][sec][bin]->SetZ( temp_max );
					}

					delete temp;

				}
				gThetaPhi[i][sec][bin] = new TGraph( nPoints, binCenter, binVal );
				for( int fun = 0; fun < 1; fun++ ){
					gThetaPhi[i][sec][bin]->Fit( Form( "fThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[i].Data() ));
				}
				
				//Draw fits
				TCanvas c(Form("hThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[i].Data()));
				hThetaPhi[i][sec][bin]->Draw("COLZ");
				gThetaPhi[i][sec][bin]->SetName(Form( "fThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[i].Data()) );
				fThetaPhi[i][sec][bin]->Write();
				fThetaPhi[i][sec][bin]->SetLineColor(kRed);
				fThetaPhi[i][sec][bin]->Draw("SAME");
				TLine * upper = new TLine( fitBounds[i][sec][bin]->X(), 0, fitBounds[i][sec][bin]->X(), 40 );
				TLine * lower = new TLine( fitBounds[i][sec][bin]->Y(), 0, fitBounds[i][sec][bin]->Y(), 40 );
				TLine * maximum = new TLine( -360, fitBounds[i][sec][bin]->Z(), 360, fitBounds[i][sec][bin]->Z() );
				upper->SetLineColor(kRed);
				lower->SetLineColor(kRed);
				maximum->SetLineColor(kRed);
				upper->Draw("SAME");
				lower->Draw("SAME");
				maximum->Draw("SAME");
				c.Write();	
				fitBounds[i][sec][bin]->Write( Form("fitBounds_sec_%i_bin_%i_%s", sec, bin, particle[i].Data() ) );
			}
		}
	}	



	/*
	for( int sec = 0; sec < 6; sec++ ){
		for( int idx = 0; idx < 2; idx++ ){
			fitMin[sec][idx] = new TF1(Form("min_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);
			fitMin[sec][idx]->SetParameters(6, 7); 
			fitMin[sec][idx]->SetParNames("a", "b"); 

			fitMax[sec][idx] = new TF1(Form("max_%i_", sec)+label[idx], "[0]+[1]/x", 1.25, 4.75);;
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
	*/

	outFile->Close();
}

double getThetaPct( double pct, TH1D * h){
	double cutoff = h->Integral()*pct;
	
	double cumCount = 0;

	for ( int i = 1; i <= h->GetNbinsX(); i++ ){
		cumCount += h->GetBinContent(i);
		
		if( cumCount >= cutoff ){
			return h->GetBinCenter(i);
		}
	}

	return h->GetBinCenter( h->GetNbinsX() );
}



