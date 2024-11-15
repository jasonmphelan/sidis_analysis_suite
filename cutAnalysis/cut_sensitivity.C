#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCut.h"
#include "TPad.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

const int num_trial = 10;
const double deg_to_rad = TMath::Pi()/180.;
const double rad_to_deg = 180./TMath::Pi();

const TLorentzVector p_rest(0., 0., 0., 0.938);

const int bins_Z = 8;
const double Z_min = 0.3;
const double Z_max = 0.8;

const double pips_parameters[6][2] = {{5.83055147, 1.16171268},
                       {5.78469108, 0.80349681},
                       {5.844136,   0.53329847},
                       {5.81600646, 0.62782215},
                       {5.75400247, 0.88392704},
                       {5.4041441,  2.12325929}};

const double pims_parameters[6][2] = {{ 6.72342364, 14.60912754},
                       { 6.85593523, 14.2618346 },
                       { 6.66919406, 14.64725313},
                       { 6.77248325, 14.33494692},
                       { 6.63926263, 14.4344934 },
                       { 6.85481318, 14.3687773 }};


std::vector<double> getCutList(double parameter, double width, int len);
double computeMean(double list[], int length);
double computeStdDev(double list[], double mean, int length);
double * computeStdDev_t(double list[], double mean);

void cut_sensitivity(){

	//Get Cut Lists	
	std::vector<double> cut_Mx_max = getCutList(5., 0.05, num_trial); 
	std::vector<double> cut_Mx_min = getCutList(1.7, 0.05, num_trial); 

	std::vector<double> cut_W_min = getCutList(2.5, 0.05, num_trial); 

	//std::vector<double> cut_Q2_min = getCutList(2., 0.05, num_trial); 
	
	std::vector<double> cut_y_max = getCutList(.75, 0.05, num_trial); 
	
	std::vector<double> cut_P_e_min = getCutList(3., 0.05, num_trial); 
	std::vector<double> cut_P_e_max = getCutList(10.6, 0.05, num_trial); 
	
	std::vector<double> cut_P_pi_min = getCutList(1.25, 0.05, num_trial); 
	std::vector<double> cut_P_pi_max = getCutList(5., 0.05, num_trial); 
	
	std::vector<double> cut_theta_e_min = getCutList(5., 0.05, num_trial); 
	std::vector<double> cut_theta_e_max = getCutList(35., 0.05, num_trial); 
	
	std::vector<double> cut_theta_pi_min = getCutList(5., 0.05, num_trial); 
	std::vector<double> cut_theta_pi_max = getCutList(35., 0.05, num_trial); 

	TFile * outFile = new TFile("cut_sensitivity_hist.root", "RECREATE");
	

	TFile * inFile = new TFile("/volatile/clas12/users/jphelan/SIDIS/data/data_no_cuts_e_piplus.root");
	TTreeReader r_plus("piplus", inFile);
	
	TTreeReaderValue<int> Npi_ptr_p(r_plus, "Npi");


        TTreeReaderValue<TLorentzVector> e_p(r_plus, "e");
        TTreeReaderValue<TLorentzVector> q_p(r_plus, "q");

        TTreeReaderValue<double> Q2_ptr_p(r_plus, "Q2");
        TTreeReaderValue<double> xB_ptr_p(r_plus, "xB");
        TTreeReaderValue<double> W_ptr_p(r_plus, "W");
        TTreeReaderValue<double> y_ptr_p(r_plus, "y");

        TTreeReaderArray<TLorentzVector> pi_p(r_plus, "pi");
        TTreeReaderArray<double> Z_p(r_plus, "Z");
        TTreeReaderArray<int> charge(r_plus, "charge");
        TTreeReaderArray<double> sector(r_plus, "pi_DC_sector");

	TTreeReaderArray<bool> passedSelectionCuts(r_plus, "PassedSelectionCuts");
	
	double pi_bins_trial[2][bins_Z][num_trial] = {0};

	int event_count = 0;

	while(r_plus.Next()){

		if(event_count%10000 == 0){cout<<"Events Analyzed: "<<event_count<<endl;}
                event_count++;

                int Npi = *Npi_ptr_p;
                if(Npi <1){continue;}

                double Q2 = *Q2_ptr_p;
                double W = *W_ptr_p;
                double xB = * xB_ptr_p;
                double y = *y_ptr_p;

		for( int trial = 0; trial < num_trial; trial++){
			
			//Apply electron cuts
			if( W < cut_W_min[trial] && W > 100. ) { continue; }
			//if( Q2 < cut_Q2_min[trial] && Q2 > 10.0) { continue; }
			if( Q2 < 2.0 && Q2 > 8.0) { continue; }
			if( xB < .1 && xB > 0.6) { continue; }
			if( y > cut_y_max[trial] ) { continue; }
			if( e_p->Vect().Mag() < cut_P_e_min[trial] || e_p->Vect().Mag() > cut_P_e_max[trial] ) { continue; }
			if( e_p->Theta()*rad_to_deg < cut_theta_e_min[trial] || e_p->Theta()*rad_to_deg > cut_theta_e_max[trial] ){ continue; }

			for( int i = 0; i < Npi; i++ ){
				if( i > 4 ){ continue; }
				if( !passedSelectionCuts[i] ){ continue; }				
	
				double M_x = (*q_p - pi_p[i] + p_rest).Mag();

				int sector_i = (int)sector[i];
				if( sector_i <1 ){continue;}
	
				double acc_map_pip = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/pi_p[i].Vect().Mag();
                        	double acc_map_pim = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/pi_p[i].Vect().Mag();

		
				int chargeIdx = (int)( charge[i]<0 );


				if ( M_x < cut_Mx_min[trial] || M_x > cut_Mx_max[trial] ) { continue; }
				if ( pi_p[i].Vect().Mag() < cut_P_pi_min[trial] || pi_p[i].Vect().Mag() > cut_P_pi_max[trial] ) { continue; }
				if ( Z_p[i] < Z_min  ||  Z_p[i] > Z_max ) { continue; }
				if ( pi_p[i].Theta()*rad_to_deg < cut_theta_pi_min[trial] || pi_p[i].Theta()*rad_to_deg > cut_theta_pi_max[trial] ){ continue; }
				if ( pi_p[i].Theta()*rad_to_deg < acc_map_pip ){ continue; }
				if ( pi_p[i].Theta()*rad_to_deg < acc_map_pim ){ continue; }

				//identify Z bin
				int this_bin_Z = (int)( ( (Z_p[i] - Z_min)/(Z_max-Z_min) )*bins_Z);
				pi_bins_trial[chargeIdx][this_bin_Z][trial]++;
			}
		}
	}
	inFile->Close();
	
	outFile->cd();

	//compute ratio for each trial
	
	double ratio[bins_Z][num_trial];
	double ratio_max[bins_Z];
	double ratio_min[bins_Z];
	double sigma_max[bins_Z];
	double sigma_min[bins_Z];
	double * sigma;
	
	cout<<"COMPUTE RATIO///////////////\n";
	for( int i = 0; i < num_trial; i++ ){
		for( int j = 0; j < bins_Z; j++ ){
			double temp_ratio = pi_bins_trial[0][j][i]/pi_bins_trial[1][j][i];
			ratio[j][i] = (4.-temp_ratio)/(4*temp_ratio-1.);
			//cout<<ratio[i][j]<<endl;
		}
	}

	TCanvas * c = new TCanvas("c", "c");

	TPad * pad[10];

	for( int bin = 0; bin < 10; bin++ ){
		if( bin < 5 ){
			pad[bin] = new TPad(Form("p_%i", bin), Form("p_%i", bin), .2*bin, .5, .2*(bin+1) , 1.);
			//pad[bin]->SetBottomMargin(0.01);
		}
		else{
			pad[bin] = new TPad(Form("p_%i", bin), Form("p_%i", bin), .2*(bin - 5), 0., .2*(bin - 4) , .5);
			//pad[bin]->SetTopMargin(0.025);
			pad[bin]->SetBottomMargin(0.25);
		}
		
		if( !( bin == 0 || bin == 5 ) ){
			pad[bin]->SetLeftMargin(0.2);
		}
		/*
		else{
			pad[bin]->SetLeftMargin(.1);
		}	
		if( !( bin == 4 || bin == 9 ) ){
			pad[bin]->SetRightMargin(0.025);
		}
		*/
		pad[bin]->Draw();
		c->cd();

	}

	TH1F * hRatio_Z[bins_Z];
	double x_min[10] = { 0.4, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1 };
	double x_max[10] = { 0.7, 0.7, 0.6, 0.6, 0.8, 0.6, 0.6, 0.6, 0.8, 0.8 };
	double maximum = 0;
	for( int i = 0; i < bins_Z; i++ ){

		hRatio_Z[i] = new TH1F(Form("hZ_%i", i), "; r(z); Number of Trials", 700, .1, .8);
		hRatio_Z[i]->SetStats(0);
		for( int j = 0; j < num_trial; j++ ){
			hRatio_Z[i]->Fill(ratio[i][j]);
		}
		hRatio_Z[i]->Write();
		double max_temp = hRatio_Z[i]->GetMaximum();
		if( max_temp > maximum ){
			maximum = max_temp;
		}
		
		hRatio_Z[i]->GetXaxis()->SetLabelSize(.055);
		hRatio_Z[i]->GetXaxis()->SetTitleSize(.075);
		
		//if( ( i == 0 || i == 5 ) ){
		hRatio_Z[i]->GetYaxis()->SetLabelSize(.055);
		hRatio_Z[i]->GetYaxis()->SetTitleSize(.075);

		//}
			
		
	}
	for( int i = 0; i < bins_Z ; i++ ){
		hRatio_Z[i]->GetYaxis()->SetRangeUser(0, 1.3*maximum);
		pad[i]->cd();
		hRatio_Z[i]->Fit("gaus");
		hRatio_Z[i]->Draw();
	}
	
	
	c->Write();
	c->SaveAs("ratio_freq.pdf");

	
	TH1F * ratio_plot_max = new TH1F("hResult_max", "; z; r(z)", 10, .3, .8);
	TH1F * ratio_plot_min = new TH1F("hResult_min", "; z; r(z)", 10, .3, .8);
	for( int i = 0; i < bins_Z; i++ ){		
		TF1 * f = hRatio_Z[i]->GetFunction("gaus");
		double mean = f->GetParameter(1);
		double sigma = f->GetParameter(2);
		ratio_plot_max->SetBinContent(i+1, mean+sigma);
		ratio_plot_min->SetBinContent(i + 1, mean-sigma);
		
	}

	TF1 * hFF = new TF1("FF", "(1-x)/(1-x+x/0.46)", Z_min, Z_max);

        TCanvas * c1 = new TCanvas("c1","c1");
        c1->Divide(1,2);
	c1->cd(1);
	ratio_plot_max->SetLineColor(kRed);
      	ratio_plot_max->GetYaxis()->SetRangeUser(0.,.75);
	ratio_plot_max->Draw();
	ratio_plot_max->Write(); 
        ratio_plot_min->SetLineColor(kBlue);
        ratio_plot_min->Draw("SAME");
	ratio_plot_max->Write();
	//hFF->SetLineColor(kViolet);
	//hFF->Draw("SAME");
	hFF->SetLineColor(kBlack);	
	hFF->Draw("SAME");

	TH1F * diff = (TH1F *)ratio_plot_max->Clone();
	diff->Add(ratio_plot_min, -1);

	c1->cd(2);
	diff->GetYaxis()->SetRangeUser(0, .05);
	ratio_plot_max->Add(ratio_plot_min);
	ratio_plot_max->Scale(1./2.);
	diff->Divide(ratio_plot_max);
	diff->GetYaxis()->SetTitle("delta R/R");
	diff->Draw();

	c1->Write();
	c1->SaveAs("cut_sensitivity.pdf");			
/*
	TH1F * hSigma = new TH1F("hSigma", "hSigma;Number of trials;#sigma", num_trial, 0, num_trial);
	
	for( int i = 0; i < num_trial; i++){
		hSigma->SetBinContent(i, *(sigma + i));
	}
	
	TCanvas * c2 = new TCanvas("c2","c2");
	hSigma->Draw();
	
	c2->SaveAs("sigma.png");
*/

	outFile->Close();	
	
}

std::vector<double> getCutList(double parameter, double width, int len){

	std::vector<double> cutList;
	std::default_random_engine generator;
	double std = parameter*width;
	std::normal_distribution<double> distribution(parameter, std);

	for( int i = 0 ; i < len ; i++ ){
		double cut_value = distribution(generator);
		cutList.push_back(cut_value);
	}

	return cutList;
}

double computeMean(double list[], int length){
	double sum_temp = 0;
	for ( int i = 0; i < length; i++ ) {
		if(TMath::IsNaN(list[i])){continue;};
		sum_temp = sum_temp + list[i];
	}
	
	double mean = sum_temp/(double)length;
	
	return mean;
}

double computeStdDev(double list[], double mean, int length){
	double num_temp = 0;

	for (int  i = 0; i < length; i ++){
		if(TMath::IsNaN(list[i])){continue;};
		num_temp = num_temp + pow(mean - list[i] , 2);
	}
	
	double var = ((double)num_temp)/((double)length);
	
	double std_dev = sqrt(var);

	return std_dev;
}

double * computeStdDev_t(double list[], double mean){
	double num_temp = 0;
	static double std_dev[num_trial];

	for (int  i = 0; i < num_trial; i ++){
		if(TMath::IsNaN(list[i])){continue;};
		num_temp = num_temp + pow(mean - list[i] , 2);
		double var = ((double)num_temp)/((double)(i+1));

		std_dev[i] = sqrt(var);
	}

	return std_dev;
}

