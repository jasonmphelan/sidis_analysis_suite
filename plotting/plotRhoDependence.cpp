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

void makeComparisonCanvas(TH1F * h1, TH1F * h2, TString varTit, TString tit_1, TString tit_2, TString outfileName){

	//h1->Sumw2();
	//h2->Sumw2();

	double fontSize = 30;
	double labelSize = 25;
	double titleSize = 30;
	int fontStyle = 43;

	h1->Scale(1./(double)h1->Integral());
	h2->Scale(1./ (double)h2->Integral());
	
	
	double temp_h1 =h1->GetMaximum();
	double temp_h2 =h2->GetMaximum();
	double maximum;
	if(temp_h1 >= temp_h2){ maximum = temp_h1; }
	else{ maximum = temp_h2; }

	TString temp_x = varTit;

	h1->GetYaxis()->SetRangeUser(0, 1.3*maximum);
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("");
	h2->GetYaxis()->SetRangeUser(0, 1.3*maximum);

	h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kAzure);
	h1->SetMarkerColor(kAzure);
	h2->SetLineColor(kRed);
	h1->SetStats(0);

	

	TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
	TPad * upper = new TPad("p1", "p1", 0, .3, 1, 1);
	upper->SetBottomMargin(0.017);
	upper->SetLeftMargin(0.15);

	upper->Draw();
	upper->cd();
	h1->Draw("E");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(0, "x");
	h1->SetLabelSize(labelSize, "y");

	h1->SetLineWidth(2);
	h2->SetLineWidth(1);

	h1->GetYaxis()->SetTitleSize(titleSize);
	h2->Draw("hist E SAME");
	
	c1->cd();
	
	TPad * lower = new TPad("p2", "p2", 0, 0.05, 1, .3);
	lower->SetTopMargin(0.017);
	lower->SetBottomMargin(0.4);
	lower->SetLeftMargin(0.15);
	lower->Draw();
	lower->cd();

	//auto ratio_Z =new TRatioPlot(h1, h2);
	//ratio_Z->SetH1DrawOpt("E1");	
	//ratio_Z->SetH2DrawOpt("E1");	
	//ratio_Z->Draw();

	TH1F * ratio = (TH1F *)h1->Clone();
	
	ratio->Divide(h2);
	ratio->Sumw2();
	ratio->SetStats(0);
	ratio->SetTitle("");

	ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	ratio->SetMarkerStyle(8);
	ratio->SetLineColor(kBlack);
	ratio->SetMarkerColor(kBlack);
	ratio->GetYaxis()->SetTitle("Ratio");
	ratio->GetXaxis()->SetTitle(temp_x);
	
	ratio->GetXaxis()->SetTitleFont(fontStyle);
	ratio->GetXaxis()->SetTitleSize(titleSize);
	ratio->GetXaxis()->SetTitleOffset(0.9);
	ratio->GetYaxis()->SetTitleFont(fontStyle);
	ratio->GetYaxis()->SetTitleSize(titleSize);
	
	ratio->GetXaxis()->SetLabelSize(labelSize);
	ratio->GetYaxis()->SetLabelSize(labelSize);
	ratio->GetYaxis()->SetNdivisions(505);


	ratio->Draw();

	TF1 * line_1 = new TF1("line_1", "1", h1->GetBinLowEdge(1), h1->GetBinLowEdge(h1->GetNbinsX() + 1) );
	line_1->SetLineColor(kBlack);
	line_1->SetLineStyle(2);
	line_1->SetLineWidth(1);

	upper->cd();
	TLegend * legend = new TLegend(0.67,0.71, .89, .89);
	legend->SetHeader("Legend", "C");
	legend->AddEntry(h2, tit_2, "l");
	legend->AddEntry(h1, tit_1, "l");
	legend->Draw();

	c1->SaveAs(outfileName);

	delete c1;
	delete ratio;
	delete h1;
	delete h2;
//	delete ratio_Z;
}

void formatHist(TH1F * h){

	double fontSize = 30;
	double labelSize = 25;
	double titleSize = 30;
	int fontStyle = 43;

	h1->Scale(1./(double)h1->Integral());
	
	
	h1->GetYaxis()->SetRangeUser(0, 1.3*maximum);
	h1->SetTitle("");
	h1->GetXaxis()->SetTitle("");

	h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kAzure);
	h1->SetMarkerColor(kAzure);
	h1->SetStats(0);

	

	h1->Draw("E");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(0, "x");
	h1->SetLabelSize(labelSize, "y");

	h1->SetLineWidth(2);

	h1->GetYaxis()->SetTitleSize(titleSize);
	
}

void plotRhoDependence(){
	TFile * file = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/rho_dependence_10.4.root");

	
	for( int x = 1; x <= 10; x++ ){
		TCanvas * q2_dep = new TCanvas("c1", "c1", 1600, 1000);
		TLegend * legend = new TLegend(0.67,0.71, .89, .89);
		legend->SetHeader("Legend", "C");
		
		for( int q2 = 1; q2 <= 12; q2++ ){

			TH1F * z_pip = (TH1F *)file->Get(Form("hZ_pip_%i_%i", q2, x));	       
			TH1F * z_r_pip = (TH1F *)file->Get(Form("hZ_r_pip_%i_%i", q2, x));	       
			TH1F * z_pim = (TH1F *)file->Get(Form("hZ_pim_%i_%i", q2, x));	       
			TH1F * z_r_pim = (TH1F *)file->Get(Form("hZ_r_pim_%i_%i", q2, x));	     

			makeComparisonCanvas(z_r_pip, z_r_pim, "z", "#pi^{+}", "#pi^{-}", Form("chargeDep_%i_%i.pdf", q2, xB) ){
			makeComparisonCanvas(z_pip, z_r_pip, "z", "SIDIS", "Diffractive #rho", Form("rhoEff_pip_%i_%i.pdf", q2, xB) ){
			makeComparisonCanvas(z_pim, z_r_pim, "z", "SIDIS", "Diffractive #rho", Form("rhoEff_pim_%i_%i.pdf", q2, xB) ){
	       			
			z_r_pip->Divide(z_pip);
			z_r_pim->Divide(z_pim);
			formatHist( z_r_pip );
			z_r_pip->SetTitle( Form( "%.2f < x_{B} < %.2f", .15 + (x-1)-
			formatHist( z_r_pim );



