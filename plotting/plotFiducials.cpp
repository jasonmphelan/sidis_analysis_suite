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


void format(TH1F * h1){

	double fontSize = 35;
	double labelSize = 30;
	double titleSize = 35;
	int fontStyle = 43;

	//h1->Scale(1./(double)h1->Integral());
	
	
	double maximum = h1->GetMaximum();


	h1->GetYaxis()->SetRangeUser(0, 1.3*maximum);
	h1->SetTitle("");
	//h1->GetXaxis()->SetTitle("");

	h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetStats(0);

	

	//h1->GetXaxis()->SetTitle(temp_x);
	//h1->GetYaxis()->SetTitle("Counts [a.u.]");
	
	h1->GetXaxis()->SetLabelSize(labelSize);
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "y");

	h1->SetLineWidth(2);

	h1->GetYaxis()->SetTitleSize(titleSize);


}

void plotFiducials(){
	TFile * f1 = new TFile( "../histograms/analysis_note/fiducials.root" );
	gStyle->SetPalette(kRainBow);
	TString lstring[5] = {"All", "5 < #theta < 12", "12 < #theta < 18", "18 < #theta < 24", "24 < #theta < 40"};
	for( int lay = 0; lay < 3; lay++ ){
		for( int sec = 0; sec < 6; sec++ ){//6
			TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
			TLegend * legend = new TLegend(0.67,0.71, .89, .89);
			legend->SetHeader("Legend", "C");
			for( int bin = 0; bin < 5; bin++ ){//5
				TH1F * h1 = (TH1F *)f1->Get(Form("hFid_e_w_reg_%i_%i_%i_pip", lay, bin, sec ) );
				format(h1);				
				h1->Draw("same PLC PMC");
				legend->AddEntry(h1, lstring[bin], "l");
				//delete h1;
			}
			legend->Draw();
			c1->SaveAs( Form("fiducials_reg_%i_sec_%i.pdf", lay, sec ) );
			delete c1;
		}
	}
	for( int lay = 0; lay < 3; lay++ ){
		for( int sec = 0; sec < 6; sec++ ){//6
			TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
			TLegend * legend = new TLegend(0.67,0.71, .89, .89);

			legend->SetHeader("Legend", "C");
			
			for( int bin = 0; bin < 5; bin++ ){//5
				TH1F * h1 = (TH1F *)f1->Get(Form("hFid_pi_w_reg_%i_%i_%i_pip", lay, bin, sec ) );
				format(h1);				
				h1->Draw("same PLC PMC");
				//delete h1;
				legend->AddEntry(h1, lstring[bin], "l");
			}
			legend->Draw();
			c1->SaveAs( Form("fiducials_pip_reg_%i_sec_%i.pdf", lay, sec ) );
			delete c1;
		}
	}
	for( int lay = 0; lay < 3; lay++ ){
		for( int sec = 0; sec < 6; sec++ ){//6
			TCanvas * c1 = new TCanvas("c1", "c1", 1600, 1000);
			TLegend * legend = new TLegend(0.67,0.71, .89, .89);

			legend->SetHeader("Legend", "C");
			
			for( int bin = 0; bin < 5; bin++ ){//5
				TH1F * h1 = (TH1F *)f1->Get(Form("hFid_pi_w_reg_%i_%i_%i_pim", lay, bin, sec ) );
				format(h1);				
				h1->Draw("same PLC PMC");
				//delete h1;
				legend->AddEntry(h1, lstring[bin], "l");
			}
			legend->Draw();
			c1->SaveAs( Form("fiducials_pim_reg_%i_sec_%i.pdf", lay, sec ) );
			delete c1;
		}
	}
	

}
