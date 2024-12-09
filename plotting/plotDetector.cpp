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
//#include "TIter.h"
#include "TKey.h"
#include "TClass.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

void drawCut(TString varTit, double yMax, TCanvas * c1);
void makeCanvas(TH1F * h1, TString varTit, TString outfileName);
void makeDetectorPlots();
void makeCanvas(TH2F * h1, TString xAxis, TString yAxis, TString outFileName);

void drawCut(TString varTit, double yMax, TCanvas * c1){
	c1->cd();

	if( varTit == "hEpcal" || varTit == "hEPcal" ){
		TLine * l1 = new TLine( .07, 0, .07, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
	}
	if( varTit == "hVz_e_pip" || varTit == "hVz_e_pim" ){
		TLine * l1 = new TLine( -5, 0, -5, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( -1, 0, -1, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hVz_pip" || varTit == "hVz_pim" ){
		TLine * l1 = new TLine( -5, 0, -5, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( -1, 0, -1, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
}
void drawCut(TString varTit, double xMin, double xMax, double yMin, double yMax, TCanvas * c1){
	c1->cd();

	if( varTit == "hPCAL_WV_pip" || varTit == "hPCAL_WV_pim" ){
		TLine * l1 = new TLine( 19, 19, xMax, 19 );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 19, 19, 19, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hEpcal_pip" || varTit == "hEpcal_pim" ){
		TLine * l1 = new TLine( .07, yMin, .07, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
	}
	if( varTit == "hSFcorr_pip" || varTit == "hSFcorr_pim" ){
		TLine * l1 = new TLine( 0, .2, .2, 0 );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
	}
	
	
}


void makeCanvas(TH1F * h1, TString varTit, TString outfileName){

	//h1->Sumw2();
	//h2->Sumw2();

	double fontSize = 35;
	double labelSize = 30;
	double titleSize = 35;
	int fontStyle = 43;

	//h1->Scale(1./(double)h1->Integral());
	
	
	double maximum = h1->GetMaximum();

	TString temp_x = varTit;

	h1->GetYaxis()->SetRangeUser(0, 1.3*maximum);
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("");

	h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetStats(0);

	

	TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
	
	
	h1->GetXaxis()->SetTitle(temp_x);
	h1->GetYaxis()->SetTitle("Counts [a.u.]");
	
	h1->GetXaxis()->SetLabelSize(labelSize);
	h1->Draw("E");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "y");

	h1->SetLineWidth(2);

	h1->GetYaxis()->SetTitleSize(titleSize);

	drawCut(h1->GetName(), 1.3*maximum, c1 );

	c1->SaveAs(outfileName);

	delete c1;
	delete h1;
}

void makeCanvas(TH2F * h1, TString xAxis, TString yAxis, TString outFileName){

	//TFile * inFile_1 = new TFile(inFileName);

	//TH2F * h1 = (TH2F *)inFile_1->Get(histName);
	h1->SetStats(0);
	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;	

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(xAxis);
	h1->GetYaxis()->SetTitle(yAxis);
	gStyle->SetPalette(kBlueGreenYellow);

	TCanvas * c1 = new TCanvas("c1", "c1");
	h1->Draw("COL");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "x");
	h1->SetLabelSize(labelSize, "y");

	double nXbins = h1->GetXaxis()->GetNbins();
	double nYbins = h1->GetYaxis()->GetNbins();
	double xMin = h1->GetXaxis()->GetBinLowEdge( 1 );
	double xMax = h1->GetXaxis()->GetBinUpEdge( nXbins );
	double yMin = h1->GetYaxis()->GetBinLowEdge( 1 );
	double yMax = h1->GetYaxis()->GetBinUpEdge( nYbins );

	h1->GetYaxis()->SetTitleSize(titleSize);
	drawCut(h1->GetName(), xMin, xMax, yMin, yMax, c1);
	c1->SetLogz();	
	c1->SaveAs(outFileName);

	delete c1;
}

void makeCanvas(TH2F * h1, TH2F * h2, TString xAxis, TString yAxis, TString outFileName){

	//TFile * inFile_1 = new TFile(inFileName);

	//TH2F * h1 = (TH2F *)inFile_1->Get(histName);
	h1->SetStats(0);
	h2->SetStats(0);
	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;	

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(xAxis);
	h1->GetYaxis()->SetTitle(yAxis);

	TCanvas * c1 = new TCanvas("c1", "c1");
	//h1->SetMaximum(1);
	h1->Draw("SCAT");
	gStyle->SetPalette(kBlueGreenYellow);
	h2->Draw("COL same");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "x");
	h1->SetLabelSize(labelSize, "y");

	double nXbins = h1->GetXaxis()->GetNbins();
	double nYbins = h1->GetYaxis()->GetNbins();
	double xMin = h1->GetXaxis()->GetBinLowEdge( 1 );
	double xMax = h1->GetXaxis()->GetBinUpEdge( nXbins );
	double yMin = h1->GetYaxis()->GetBinLowEdge( 1 );
	double yMax = h1->GetYaxis()->GetBinUpEdge( nYbins );

	h1->GetYaxis()->SetTitleSize(titleSize);
	drawCut(h1->GetName(), xMin, xMax, yMin, yMax, c1);
	c1->SetLogz();	
	c1->SaveAs(outFileName);

	delete c1;
}

void plotDetector(){
	TFile * f1 = new TFile( "test.root" );
	TIter keyList( f1->GetListOfKeys() );
	TKey *key;
	
	while( (key = (TKey*)keyList()) ){
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if( cl->InheritsFrom("TH2")){
			TH2F * h1 = (TH2F*)key->ReadObj();
			//TH2F * h1 = (TH2F *)f1->Get("hSF_sec_0_pip");
			makeCanvas(h1, "VAR", "VAR", (TString)h1->GetName()+".pdf");	
		}
		else{
			TH1F * h1 = (TH1F *)key->ReadObj();
			//TH1F * h1 = (TH1F *)f1->Get("hSF_sec_1_pim");
			makeCanvas(h1, "#chi", (TString)h1->GetName() + ".pdf");
		}
	}
}
