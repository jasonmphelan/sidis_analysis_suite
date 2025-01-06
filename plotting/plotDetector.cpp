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
		
		TLine * l2 = new TLine( 1, 0, 1, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hZ_pip" || varTit == "hZ_pim" ){
		TLine * l1 = new TLine( .3, 0, .3, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
	}
	if( varTit == "hY_pip" || varTit == "hY_pim" ){
		TLine * l1 = new TLine( .75, 0, .75, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
	}
	if( varTit == "hW_pip" || varTit == "hW_pim" ){
		TLine * l1 = new TLine( 2.5, 0, 2.5, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
	}
	if( varTit == "hTheta_e_pip" || varTit == "hTheta_e_pim" ){
		TLine * l1 = new TLine( 5, 0, 5, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 35, 0, 35, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hTheta_pi_pip" || varTit == "hTheta_pi_pim" ){
		TLine * l1 = new TLine( 5, 0, 5, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 35, 0, 35, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hQ2_pip" || varTit == "hQ2_pim" ){
		TLine * l1 = new TLine( 2, 0, 2, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
	}
	if( varTit == "hP_pi_pip" || varTit == "hP_pi_pim" ){
		TLine * l1 = new TLine( 1.25, 0, 1.25, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 5, 0, 5, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hP_e_pip" || varTit == "hP_e_pim" ){
		TLine * l1 = new TLine( 3, 0, 3, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 10, 0, 10, yMax );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
	}
	if( varTit == "hMx_pip" || varTit == "hMx_pim" ){
		TLine * l1 = new TLine( 1.7, 0, 1.7, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
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
	if( varTit == "hEdep_pip" || varTit == "hEdep_pim" ){
		TLine * l1 = new TLine( .07, yMin, .07, yMax );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
	}
	if( varTit == "hSF_corr_pip" || varTit == "hSF_corr_pim" ){
		TLine * l1 = new TLine( 0, .2, .2, 0 );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
	}
	if( varTit == "hChi2_pip" || varTit == "hChi2_pim" ){
		TLine * l1 = new TLine( 0, -3, xMax, -3 );
		l1->SetLineColor( kRed );
		l1->SetLineWidth(2);	
		l1->Draw("same");
		
		TLine * l2 = new TLine( 0, 3, 2.44, 3 );
		l2->SetLineColor( kRed );
		l2->SetLineWidth(2);	
		l2->Draw("same");
		
		TF1 * f1 = new TF1("f1", "[0]+[1]*exp(-x/[2]) + [3]*exp(-x/[4])", 2.44, 4.6);
		f1->SetParameters(0.00869, 14.98587, 1.18236, 1.81751, 4.86394);
		f1->SetLineColor( kRed );
		f1->SetLineWidth(2);	
		f1->Draw("same");
		
		TF1 * f2 = new TF1("f2", "[0]+[1]*exp(-x/[2]) + [3]*exp(-x/[4])", 4.6, xMax);
		f2->SetParameters(-1.14099,24.14992,1.36554,2.66876,6.80552 );
		f2->SetLineColor( kRed );
		f2->SetLineWidth(2);	
		f2->Draw("same");
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
	//h1->GetXaxis()->SetTitle("");

	h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetStats(0);

	

	TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
	
	
	//h1->GetXaxis()->SetTitle(temp_x);
	//h1->GetYaxis()->SetTitle("Counts [a.u.]");
	
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
	//h1->GetXaxis()->SetTitle(xAxis);
	//h1->GetYaxis()->SetTitle(yAxis);
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
	//h1->GetXaxis()->SetTitle(xAxis);
	//h1->GetYaxis()->SetTitle(yAxis);

	TCanvas * c1 = new TCanvas("c1", "c1");
	//h1->SetMaximum(1);
	h1->Draw("SCAT");
	gStyle->SetPalette(kBird);
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
	TFile * f1 = new TFile( "../histograms/analysis_note/detector_plots_10.6.root" );
	TIter keyList( f1->GetListOfKeys() );
	TKey *key;
	
	while( (key = (TKey*)keyList()) ){
		TClass *cl = gROOT->GetClass(key->GetClassName());
		bool madeHist = false;
		if( cl->InheritsFrom("TH2")){
			TH2F * h1 = (TH2F*)key->ReadObj();
			for( int i = 0; i < 3; i++ ){
				if( (TString)h1->GetName() == Form("hFid_pi_bef_reg_%i_pip", i) ){
					TH2F * h2 = (TH2F*)f1->Get(Form("hFid_pi_aft_reg_%i_pip", i));
					makeCanvas(h2, h1, "VAR", "VAR", "selection_plots/"+(TString)h1->GetName()+".pdf");	
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_pi_aft_reg_%i_pip", i) ){
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_pi_bef_reg_%i_pim", i) ){
					TH2F * h2 = (TH2F*)f1->Get(Form("hFid_pi_aft_reg_%i_pip", i));
					makeCanvas(h2, h1, "VAR", "VAR", "selection_plots/"+(TString)h1->GetName()+".pdf");	
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_pi_aft_reg_%i_pim", i) ){
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_e_bef_reg_%i_pip", i) ){
					TH2F * h2 = (TH2F*)f1->Get(Form("hFid_e_aft_reg_%i_pip", i));
					makeCanvas(h2, h1, "VAR", "VAR", "selection_plots/"+(TString)h1->GetName()+".pdf");	
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_e_aft_reg_%i_pip", i) ){
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_e_bef_reg_%i_pim", i) ){
					TH2F * h2 = (TH2F*)f1->Get(Form("hFid_e_aft_reg_%i_pip", i));
					makeCanvas(h2, h1, "VAR", "VAR", "selection_plots/"+(TString)h1->GetName()+".pdf");	
					madeHist = true;
				}
				if( (TString)h1->GetName() == Form("hFid_e_aft_reg_%i_pim", i) ){
					madeHist = true;
				}
			}

			if( madeHist ){continue;}

			makeCanvas(h1, "VAR", "VAR", "selection_plots/"+(TString)h1->GetName()+".pdf");	

		}
		else{
			TH1F * h1 = (TH1F *)key->ReadObj();
			//TH1F * h1 = (TH1F *)f1->Get("hSF_sec_1_pim");
			makeCanvas(h1, "#chi", "selection_plots/"+(TString)h1->GetName() + ".pdf");
		}
	}
}
