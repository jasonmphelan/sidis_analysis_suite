#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRatioPlot.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

void formatHist2D(TH2F * h1, TString xAxis, TString yAxis){

	h1->SetStats(0);
	double fontSize = 35;
	double labelSize = 30;
	double titleSize = 35;
	int fontStyle = 43;	

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(xAxis);
	h1->GetYaxis()->SetTitle(yAxis);

	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "x");
	h1->SetLabelSize(labelSize, "y");

	h1->GetYaxis()->SetTitleSize(titleSize);
	
}


void acceptanceMatching( TString inFileName, TString outFileDir ){
		
	TString chargeStr[2] = {"pip", "pim"};
	TString chargeType[2] = {"+", "-"};
	TFile * inFile = new TFile(inFileName);
	for( int i = 0; i < 6; i++ ){ //Sectors
		TF1 * max_pip = (TF1*) inFile->Get( Form( "max_%i_pip", i ) );
		TF1 * max_pim = (TF1*) inFile->Get( Form( "max_%i_pim", i ) );
		TF1 * min_pip = (TF1*) inFile->Get( Form( "min_%i_pip", i ) );
		TF1 * min_pim = (TF1*) inFile->Get( Form( "min_%i_pim", i ) );
	
		max_pip->SetLineColor( kGreen );
		min_pip->SetLineColor( kGreen );
		
		max_pip->SetLineWidth( 2 );
		min_pip->SetLineWidth( 2 );
		max_pim->SetLineWidth( 2 );
		min_pim->SetLineWidth( 2 );
	
		for( int k = 0; k < 2; k++ ){ //Charge
			TCanvas * c1 = new TCanvas( "c1", "c1", 1600, 800 );		
			TH2F * h = (TH2F*) inFile->Get( Form("hTheta_P_sec_%i_", i) + chargeStr[k] );
			
			formatHist2D( h, "p_{#pi} [GeV]", "#theta_{#pi} [deg.]" );
			h->SetTitle("(e, e'#pi"+chargeType[k]+"), Sector "+Form("%i", i+1) );
			h->Draw("COL");
			max_pip->Draw("SAME");
			min_pip->Draw("SAME");
			max_pim->Draw("SAME");
			min_pim->Draw("SAME");
			
			c1->SaveAs( outFileDir + Form("/acc_match_2d_%i_", i) + chargeStr[k] + ".pdf");
			delete c1;
			
		}
	}
}

void make2DCanvas(TString inFileName, TString outFileName, TString histName, TString xAxis, TString yAxis){

	TFile * inFile_1 = new TFile(inFileName);

	TH2F * h1 = (TH2F *)inFile_1->Get(histName);
	h1->SetStats(0);
	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;	

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(xAxis);
	h1->GetYaxis()->SetTitle(yAxis);

	TCanvas * c1 = new TCanvas("c1", "c1");

	h1->Draw("COLZ");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "x");
	h1->SetLabelSize(labelSize, "y");

	h1->GetYaxis()->SetTitleSize(titleSize);
	
	c1->SaveAs(outFileName + "/" + histName +".pdf");

	delete c1;
}
