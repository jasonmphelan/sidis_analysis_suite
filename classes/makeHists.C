#include "makeHists.h"


void makeCanvas(TString inFileName_1,TString outFileName, TString var, TString varTit){

	TFile * inFile_1 = new TFile(inFileName_1);

	TH1F * h1 = (TH1F *)inFile_1->Get("h"+var);
	h1->Sumw2();

	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;

	TString temp_x = varTit;
	TCanvas * c1 = new TCanvas("c1", "c1");

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(varTit);
	h1->GetYaxis()->SetTitle("Counts [a.u.]");

	//h1->SetMarkerStyle(kFullCircle);
	h1->SetLineColor(kAzure);
	h1->SetMarkerColor(kAzure);
	h1->SetStats(0);
	

	h1->Draw("hist");
	gPad->SetLogy();
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "xy");

	h1->SetLineWidth(1);

	h1->GetYaxis()->SetTitleSize(titleSize);
	
	c1->SaveAs(outFileName);

	delete c1;
}

void makeComparisonCanvas(TString inFileName_1, TString inFileName_2,TString outFileName, TString var, TString varTit,  TString pionType_1, TString pionType_2, TString tit_1, TString tit_2){

	TFile * inFile_1 = new TFile(inFileName_1);
	TFile * inFile_2 = new TFile(inFileName_2);

	TH1F * h1 = (TH1F *)inFile_1->Get("h"+var+"_"+pionType_1);
	TH1F * h2 = (TH1F *)inFile_2->Get("h"+var+"_"+pionType_2);
	//h1->Sumw2();
	//h2->Sumw2();
	makeComparisonCanvas(h1, h2, outFileName, var, varTit,  pionType_1, pionType_2, tit_1, tit_2){
}

void makeComparisonCanvas(TH1F* h1, TH1F* h2,TString outFileName, TString var, TString varTit, TString tit_1, TString tit_2){

	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;

	//h1->Scale(1./(double)h1->Integral());
	//h2->Scale(1./ (double)h2->Integral());
	
	
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

	

	TCanvas * c1 = new TCanvas("c1", "c1");
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
	TLegend * legend = new TLegend(0.67,0.71, .89, .89);
	legend->SetHeader("Legend", "C");
	legend->AddEntry(h2, tit_2, "l");
	legend->AddEntry(h1, tit_1, "l");
	legend->Draw();

	c1->SaveAs(outFileName);

	delete c1;
	delete ratio;
	delete h1;
	delete h2;
	delete inFile_1;
	delete inFile_2;
	delete ratio_Z;
}

void makeComparisonCanvas(TString inFileName_1, TString inFileName_2,TString outFileName, TString var, TString varTit,  TString pionType, TString tit_1, TString tit_2){

	makeComparisonCanvas(inFileName_1, inFileName_2, outFileName, var, varTit,  pionType, pionType, tit_1, tit_2);

}

void makeAllComparisonChargeCanvas(TString inFileName_1, TString inFileName_2, int nBinsQ2, int nBinsXb, TString outFileBase){

	std::vector<TString> variables_e{"W", "Xb", "Q2", "Y", "Omega", "Pt_e", "Vz_e", "P_e", "Phi_e"}; //Theta_e
	std::vector<TString> variable_titles_e{"W", "x_{B}", "Q^{2}", "y", "#omega", "p^{#perp}_{e}", "v^{z}_{e}", "p_{e}", "#phi_{e}"};
	
	std::vector<TString> variables_pi{"Z", "P_pi", "Vz_pi", "DeltaPhi", "Theta_pi", "Phi_pi", "Pt_pi", "Mx"};
	std::vector<TString> variable_titles_pi{"z", "p_{#pi}", "v_{#pi}^{z}", "#phi_{e} - #phi_{#pi}","#theta_{#pi}", "#phi_{#pi}", "p^{#perp}_{#pi}", "M_{X}"};
		
	for( int k = 2; k <= 5; k++ ){
		for( int j = 2; j <= 5; j++ ){	
			//for( int i = 0; i < 9; i++ ){
			//	makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase + variables_e[i] +Form("_pip_%i_%i.pdf", j, k), variables_e[i], variable_titles_e[i],  Form("pip_%i_%i", j, k), Form("pim_%i_%i", j, k) , "#pi^{+}", "#pi^{-}");
			//}
			for( int i = 0; i < 8; i++ ){
				makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase + variables_pi[i] +Form("_%i_%i.pdf", j, k), variables_pi[i], variable_titles_pi[i],  Form("pip_%i_%i", j, k), Form("pim_%i_%i", j, k) , "#pi^{+}", "#pi^{-}");
			}
		}
	}

}

void makeAllComparisonCanvas(TString inFileName_1, TString inFileName_2, int nBinsQ2, int nBinsXb, TString outFileBase, TString tit1, TString tit2){

	std::vector<TString> variables_e{"W", "Xb", "Q2", "Y", "Omega", "Pt_e", "Vz_e", "P_e", "Phi_e"}; //Theta_e
	std::vector<TString> variable_titles_e{"W", "x_{B}", "Q^{2}", "y", "#omega", "p^{#perp}_{e}", "v^{z}_{e}", "p_{e}", "#phi_{e}"};
	
	std::vector<TString> variables_pi{"Z", "P_pi", "Vz_pi", "DeltaPhi", "Theta_pi", "Phi_pi", "Pt_pi", "Mx"};
	std::vector<TString> variable_titles_pi{"z", "p_{#pi}", "v_{#pi}^{z}", "#phi_{e} - #phi_{#pi}","#theta_{#pi}", "#phi_{#pi}", "p^{#perp}_{#pi}", "M_{X}"};
		
	for( int k = 2; k <= 5; k++ ){
		for( int j = 2; j <= 5; j++ ){	
			for( int i = 0; i < 9; i++ ){
				makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase + variables_e[i] +Form("_pip_%i_%i.pdf", j, k), variables_e[i], variable_titles_e[i],  Form("pip_%i_%i", j, k) , tit1, tit2);
				makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase + variables_e[i] +Form("_pim_%i_%i.pdf", j, k), variables_e[i], variable_titles_e[i],  Form("pim_%i_%i", j, k) , tit1, tit2);
			}
			for( int i = 0; i < 8; i++ ){
				makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase +  variables_pi[i] +Form("_pip_%i_%i.pdf", j, k), variables_pi[i], variable_titles_pi[i],  Form("pip_%i_%i", j, k), tit1, tit2);
				makeComparisonCanvas("rootFiles/kinematics/"+inFileName_1, "rootFiles/kinematics/"+inFileName_2, outFileBase + variables_pi[i] +Form("_pim_%i_%i.pdf", j, k), variables_pi[i], variable_titles_pi[i],  Form("pim_%i_%i", j, k), tit1, tit2);
			}
		}
	}

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void make2DCanvasLog(TString inFileName, TString outFileName, TString histName, TString xAxis, TString yAxis){
	TFile * inFile_1 = new TFile(inFileName);

	TH2F * h1 = (TH2F *)inFile_1->Get(histName);

	double fontSize = 20;
	double labelSize = 15;
	double titleSize = 20;
	int fontStyle = 43;	

	h1->SetTitle("");
	h1->GetXaxis()->SetTitle(xAxis);
	h1->GetYaxis()->SetTitle(yAxis);
	h1->SetStats(0);

	TCanvas * c1 = new TCanvas("c1", "c1");

	h1->Draw("COLZ");
	h1->GetXaxis()->SetTitleFont(fontStyle);
	h1->GetXaxis()->SetTitleSize(titleSize);
	h1->GetXaxis()->SetTitleOffset(0.9);
	h1->GetYaxis()->SetTitleFont(fontStyle);
	h1->GetYaxis()->SetTitleSize(titleSize);
	gPad->SetLogz();

	h1->SetLabelFont(fontStyle, "xyz");
	h1->SetLabelSize(labelSize, "x");
	h1->SetLabelSize(labelSize, "y");

	h1->GetYaxis()->SetTitleSize(titleSize);
	
	c1->SaveAs(outFileName + "/" + histName+".pdf");

	delete c1;
}

void loopThrough2DCanvas(TString inFileName, TString outFileName){
	
	std::vector<TString> variables_pi{"Z_Mx", "Q2_Mx", "Xb_Mx", "Pt_Mx"};
	std::vector<TString> variable_titles_pi_x{"z", "Q^{2} [GeV^{2}]", "x_{B}", "P^{#perp}_{#pi}"};
	
	for( int i = 0; i < 4; i++ ){
		for( int j = 0; j <= 12; j++ ){	
			make2DCanvas("rootFiles/kinematics/"+inFileName+".root", "PDFs/kinematics/"+outFileName+"/"+variables_pi[i]+Form("_pip_%i.pdf", j), "h"+variables_pi[i]+Form("_pip_%i", j), variable_titles_pi_x[i], "M_{X} [GeV]");
			make2DCanvas("rootFiles/kinematics/"+inFileName+".root", "PDFs/kinematics/"+outFileName+"/"+variables_pi[i]+Form("_pim_%i.pdf", j), "h"+variables_pi[i]+Form("_pim_%i", j), variable_titles_pi_x[i], "M_{X} [GeV]");
		}
	}
}

void makeAllDetectorPlots(){
	
	TString histNames[10] = {"Pos", "Min", "Pipfid", "Pipvertex", "Pimfid", "ElWV", "ElminE", "ElSamplingFraction", "Pimvertex", "Elfid"};
	TString histNames2D[10] = {"Pim_b_p_vertex", "Pip_b_p_vertex", "Pos_b_p", "Min_b_p", "El_b_p_minE", "El_b_p_WV", "El_b_p_SamplingFraction", "El_b_p_fid", "Pip_b_p_fid", "Pim_b_p_fid"};
	
	for( int i = 0; i < 10; i++ ){
		makeCanvas("rootFiles/skimFiles/det_skim_hists.root", "PDFs/pid/"+histNames[i]+".pdf", histNames[i], "#beta");
		make2DCanvasLog("rootFiles/skimFiles/det_skim_hists.root", "PDFs/pid/", "h"+histNames2D[i], "p_{#pi} [GeV]", "#beta");
	}
}
