#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TH1.h"
#include "TF1.h"
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

void makeComparisonCanvas(TString inFileName_1, TString inFileName_2,TString outFileName, TString charge, int binQ, int binX ){

	TFile * inFile_1 = new TFile(inFileName_1);
	TFile * inFile_2 = new TFile(inFileName_2);

//	TH1F * h1 = (TH1F *)inFile_1->Get("h"+var+"_"+pionType_1);
//	TH1F * h2 = (TH1F *)inFile_2->Get("h"+var+"_"+pionType_2);
	

	TF1 * f1 = (TF1 *)inFile_1->Get(Form("fit%s_%i_%i", charge.Data(), binQ, binX));
	TF1 * f2 = (TF1 *)inFile_2->Get(Form("fit%s_%i_%i", charge.Data(), binQ, binX));
	//h1->Sumw2();
	//h2->Sumw2();
	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;

	//h1->Scale(1./(double)h1->Integral());
	//h2->Scale(1./ (double)h2->Integral());
	
	
	//double temp_h1 =h1->GetMaximum();
	//double temp_h2 =h2->GetMaximum();
	double maximum = 5;
	//if(temp_h1 >= temp_h2){ maximum = temp_h1; }
	//else{ maximum = temp_h2; }

	TString temp_x = "z";

	f1->GetYaxis()->SetRangeUser(0, maximum);
	f1->SetTitle("");
	f1->GetXaxis()->SetTitle("z");
	f2->GetYaxis()->SetRangeUser(0, maximum);

	//h1->SetMarkerStyle(kFullCircle);
	//h1->SetLineColor(kAzure);
	f1->SetLineColor(kAzure);
	//h1->SetMarkerColor(kAzure);
	//h2->SetLineColor(kRed);
	f2->SetLineColor(kRed);
	//h1->SetStats(0);

	
	cout<<"PREP FOR CANVAS\n";
	TCanvas * c1 = new TCanvas("c1", "c1");
	//TPad * upper = new TPad("p1", "p1", 0, .3, 1, 1);
	//upper->SetBottomMargin(0.017);
	//upper->SetLeftMargin(0.15);

	//upper->Draw();
	//upper->cd();
	//h1->Draw("E");
	f1->GetXaxis()->SetTitleFont(fontStyle);
	f1->GetXaxis()->SetTitleSize(titleSize);
	f1->GetXaxis()->SetTitleOffset(0.9);
	f1->GetYaxis()->SetTitleFont(fontStyle);
	f1->GetYaxis()->SetTitleSize(titleSize);

	f1->GetXaxis()->SetLabelFont(fontStyle);
	//f1->GetXaxis()->SetLabelSize(0, "x");
	f1->GetXaxis()->SetLabelSize(labelSize);

	//h1->SetLineWidth(1);
	//h2->SetLineWidth(1);
	
	f1->SetLineWidth(2);
	f2->SetLineWidth(2);

	f1->GetYaxis()->SetTitleSize(titleSize);
	//h2->Draw("E SAME");
	f1->Draw("");

	f2->Draw("SAME");
	//c1->cd();
/*	
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

	ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
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


	upper->cd();
*/
//	TLegend * legend = new TLegend(0.67,0.71, .89, .89);
//	legend->SetHeader("Legend", "C");
//	legend->AddEntry(f1, "10.2 [GeV]", "l");
//	legend->AddEntry(f2, "10.2 [GeV]", "l");
//	legend->Draw();

	c1->SaveAs(outFileName);

	delete c1;
	//delete ratio;
	//delete h1;
	//delete h2;
	delete f1;
	delete f2;
	delete inFile_1;
	delete inFile_2;
//	delete ratio_Z;
}

void plotFits( TString infile_1, TString infile_2 ){

	for( int x = 0; x < 14; x++ ){
		for( int q = 0; q < 12; q++ ){
			std::cout<<"Plotting : x= "<<x<<" and  q = "<<q<<std::endl;
			makeComparisonCanvas(infile_1, infile_2, Form("correctionFits/fit_pip_%i_%i.pdf", q, x), "Pip", q, x );
			makeComparisonCanvas(infile_1, infile_2, Form("correctionFits/fit_pim_%i_%i.pdf", q, x), "Pim", q, x );
		}
	}
}
