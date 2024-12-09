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

void plotCharge( TString inFileName, TString outFilePath ){

	
	std::vector<TString> variables_e{"W", "Xb", "Q2", "Y", "Omega", "Vz_e", "P_e", "Phi_e", "Theta_e"}; //Theta_e
	std::vector<TString> variable_titles_e{"W [GeV]", "x_{B}", "Q^{2} [GeV^{2}]", "y", "#omega [GeV]", "v^{z}_{e} [cm]", "p_{e} [GeV.]", "#phi_{e}", "#theta_{e} [deg.]"};
	
	std::vector<TString> variables_pi{"Z", "P_pi", "Vz_pi", "Theta_pi", "Phi_pi", "Pt_pi", "Mx"};
	std::vector<TString> variable_titles_pi{"z", "p_{#pi} [GeV]", "v_{#pi}^{z} - v_{e}^{z} [cm]","#theta_{#pi} [deg.]", "#phi_{#pi}", "p^{#perp}_{#pi} [GeV]", "M_{X} [GeV]"};
	
	TFile * inFile = new TFile( inFileName );

	for( int i = 0; i < variables_e.end() - variables_e.begin(); i ++ ){
		TH1F * hPip = (TH1F *)inFile->Get("h" + variables_e[i] + "_pip_0_0");
		TH1F * hPim = (TH1F *)inFile->Get("h" + variables_e[i] + "_pim_0_0");
		
		makeComparisonCanvas(hPip, hPim, variable_titles_e[i], "(e, e'#pi+)", "(e, e'#pi-)", outFilePath + "/charge_"+variables_e[i]+".pdf");
			
	}
	for( int i = 0; i < variables_pi.end() - variables_pi.begin(); i ++ ){
		TH1F * hPip = (TH1F *)inFile->Get("h" + variables_pi[i] + "_pip_0_0");
		TH1F * hPim = (TH1F *)inFile->Get("h" + variables_pi[i] + "_pim_0_0");
		
		makeComparisonCanvas(hPip, hPim, variable_titles_pi[i], "(e, e'#pi+)", "(e, e'#pi-)", outFilePath + "/charge_"+variables_pi[i]+".pdf");
			
	}
}	

void plotMC( TString inFileName_dat, TString inFileName_mc, TString outFilePath ){

	
	std::vector<TString> variables_e{"W", "xB", "Q2", "y", "omega", "Vz_e", "P_e", "Phi_e", "Theta_e"}; //Theta_e
	std::vector<TString> variable_titles_e{"W [GeV]", "x_{B}", "Q^{2} [GeV^{2}]", "y", "#omega [GeV]", "v^{z}_{e} [cm]", "p_{e} [GeV.]", "#phi_{e}", "#theta_{e} [deg.]"};
	
	std::vector<TString> variables_pi{"Z", "P_pi", "Vz_pi", "Theta_pi", "Phi_pi", "Pt_pi", "Mx"};
	std::vector<TString> variable_titles_pi{"z", "p_{#pi} [GeV]", "v_{#pi}^{z} - v_{e}^{z} [cm]","#theta_{#pi} [deg.]", "#phi_{#pi}", "p^{#perp}_{#pi} [GeV]", "M_{X} [GeV]"};
	
	TFile * inFile = new TFile( inFileName_dat );
	TFile * inFile_mc = new TFile( inFileName_mc );

	for( int j = 0; j < 2; j++ ){
		for( int i = 0; i < variables_e.end() - variables_e.begin(); i ++ ){
			TH1F * hDat = (TH1F *)inFile->Get("h");
			TH1F * hMC = (TH1F *)inFile_mc->Get("h");
		
			makeComparisonCanvas(hDat, hMC, variable_titles_e[i], "Data", "GEMC", outFilePath + "/charge_"+variables_e[i]+".pdf");
			
		}
		for( int i = 0; i < variables_pi.end() - variables_pi.begin(); i ++ ){
			
			TH1F * hDat = (TH1F *)inFile->Get("h");
			TH1F * hMC = (TH1F *)inFile_mc->Get("h");
		
			makeComparisonCanvas(hDat, hMC, variable_titles_pi[i], "Data", "GEMC", outFilePath + "/data_sim_"+variables_pi[i]+".pdf");
			
		}
	}
}	

