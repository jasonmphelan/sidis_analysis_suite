
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
#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<clashit>+;
#pragma link C++ class vector<clashit>+;
#endif
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

	//TF1 * hFit[3][5][6][2];	 //3 for particles, 5 for p bins, 6 sectors, 2 for max/min

	TH2F * hThetaPhi[3][6][10];
	TGraph * gThetaPhi[3][6][10];
	TF1 * fThetaPhi[3][6][10][2];

	//TString momentum_functions[2] = {"[0] + exp([1]*x+[2])", "[0] + exp(-1*[1]*x-[2])"};
	TString momentum_functions[2] = {"[0]*(x - [1])*( x-[1] ) + [2]", "[0]+(x - [1])*(x-[1])+[2]"};
	TString boundary[2] = {"fit","lower"};
	TString particle[3] = {"e","pip", "pim"};
	
	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 0; bin < 10; bin++ ){
			for( int p = 0; p < 3; p++ ){
				if(p > 0 && bin >4){ continue; }
				hThetaPhi[p][sec][bin] = (TH2F *)inFile->Get( Form("hThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[p].Data()));
				double histMean = hThetaPhi[p][sec][bin]->GetMean(1);
				double histStd = hThetaPhi[p][sec][bin]->GetRMS(1);
				double histMeanY = hThetaPhi[p][sec][bin]->GetMean(2);
				double histStdY = hThetaPhi[p][sec][bin]->GetRMS(2);
				double histMin = histMeanY - 3*histStdY;
				
				//double a = histMin - histMean*histMean*histStd;
				//double b = -10./histMean;
				//double c = log(2*histMean*histMean*histStd);

				cout<<"Min : "<<histMin<<" Width : "<<histStd<< " Mean : "<<histMean<<std::endl;

				for( int fun = 0; fun < 2; fun++ ){
					fThetaPhi[p][sec][bin][fun] = new TF1(Form("fThetaPhi_sec_%i_bin_%i_%s_%s", sec, bin, particle[p].Data(), boundary[fun].Data()), 
								momentum_functions[fun], histMean - 3*histStd, histMean + (3*histStd));
					fThetaPhi[p][sec][bin][fun]->SetParameters( 1./histStd, histMean, histMin );
			
				}
			}

		}
			
	}
	cout<<"DECLARED RELEVANT OBJECTS\n";
		
	//Get Cut Value (cut out 10% percentile)
		

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

				for( int phiBin = 0; phiBin < nPhiBins; phiBin++ ){
					TH1D * temp = new TH1D( "temp", "", hThetaPhi[i][sec][bin]->GetYaxis()->GetNbins(), 0, 40 );
					temp = (TH1D *)hThetaPhi[i][sec][bin]->ProjectionY( "temp", phiBin ,  phiBin   );
					if( temp->Integral() < 500 ){
						delete temp;
						continue;
					}
					binCenter[nPoints] =  hThetaPhi[i][sec][bin]->GetXaxis()->GetBinCenter( phiBin )  ; 
					binVal[nPoints] = getThetaPct( .05, temp ) ;
					nPoints++;
					delete temp;

				}
				gThetaPhi[i][sec][bin] = new TGraph( nPoints, binCenter, binVal );
				for( int fun = 0; fun < 1; fun++ ){
					gThetaPhi[i][sec][bin]->Fit( Form( "fThetaPhi_sec_%i_bin_%i_%s_%s", sec, bin, particle[i].Data(), boundary[fun].Data() ));
				}
				TCanvas c(Form("hThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[i].Data()));
				gThetaPhi[i][sec][bin]->Draw();
				gThetaPhi[i][sec][bin]->SetName(Form( "fThetaPhi_sec_%i_bin_%i_%s", sec, bin, particle[i].Data()) );
				gThetaPhi[i][sec][bin]->Write();
				fThetaPhi[i][sec][bin][0]->SetLineColor(kRed);
				fThetaPhi[i][sec][bin][0]->Write();
				//fThetaPhi[i][sec][bin][1]->SetLineColor(kRed);
				fThetaPhi[i][sec][bin][0]->Draw("SAME");
				//fThetaPhi[i][sec][bin][1]->Draw("SAME");
				//c.Write();	
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
		
		if( cumCount > cutoff ){
			return h->GetBinCenter(i);
		}
	}

	return h->GetBinCenter( h->GetNbinsX() );
}

