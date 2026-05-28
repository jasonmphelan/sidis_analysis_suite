#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "analyzer.h"
#include "DCfid_SIDIS.h"
#define HIST_PATH _HIST

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut );


int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		cerr << "[Type (0 - pi2k, 1 - k2pi)] [Theta Cut (optional)] [Mx_2pi^2 cut (optional, default 1.25)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";


	int p_bins = 30;

    TString out_name = argv[1];
		
	int nPhotons = atoi(argv[2]);
	double prob = atof(argv[3]);
	int data_sim = atoi(argv[4]);

    TFile * outFile = new TFile(out_name, "RECREATE");//(TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
	
	// Declare histograms

	cout<<"Creating Histograms\n";

	// 4 momentum histograms — one per species
	// Index: 0=pip misid, 1=pim misid, 2=kp, 3=km
	TString species[6] = {"pip", "pim", "kp", "km", "prot", "aprot"};
	TH1F * hP_eb[6];
	TH1F * hP_true[6];
	TH1F * hP_opp[6];


	TH2F * hTheta_p_eb[6];
	TH2F * hTheta_p_true[6];
	TH2F * hNRICH_p[6];
	TH2F * hProbRICH_p[6];
	TH2F * hNRICH_p_eb[6];
	TH2F * hProbRICH_p_eb[6];

	for( int i = 0; i < 6; i++ ){
		hP_eb[i] = new TH1F(
			"hP_eb_" + species[i],
			"p (" + species[i] + ");p (GeV/c);Counts",
			p_bins, 1.25, 5.0 );
		hP_true[i] = new TH1F(
			"hP_true_" + species[i],
			"p (" + species[i] + ");p (GeV/c);Counts",
			p_bins, 1.25, 5.0 );
		hTheta_p_eb[i] = new TH2F(
			"hTheta_p_eb_" + species[i],
			"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
			100, 1.25, 5.0, 90, 5.0, 35.0 );
		hTheta_p_true[i] = new TH2F(
			"hTheta_p_true_" + species[i],
			"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
			100, 1.25, 5.0, 90, 5.0, 35.0 );
		hP_opp[i] = new TH1F(
				"hP_opp_" + species[i],
				"p (" + species[i] + ");p (GeV/c);Counts",
				p_bins, 1.25, 5.0 );
		hNRICH_p[i] = new TH2F(
				"hNRICH_p_" + species[i],
				"n_{RICH} vs p (" + species[i] + ");p (GeV/c);n_{RICH}",
				100, 1.25, 5.0, 50, 0, 50 );
		hProbRICH_p[i] = new TH2F(
				"hProbRICH_p_" + species[i],
				"prob_{RICH} vs p (" + species[i] + ");p (GeV/c);prob_{RICH}",
				100, 1.25, 5.0, 100, 0, 1.0 );
		hNRICH_p_eb[i] = new TH2F(
				"hNRICH_p_eb_" + species[i],
				"n_{RICH} vs p, RICH-matched (" + species[i] + ");p (GeV/c);n_{RICH}",
				100, 1.25, 5.0, 50, 0, 50 );
		hProbRICH_p_eb[i] = new TH2F(
				"hProbRICH_p_eb_" + species[i],
				"prob_{RICH} vs p, RICH-matched (" + species[i] + ");p (GeV/c);prob_{RICH}",
				100, 1.25, 5.0, 100, 0, 1.0 );
	}

	cout<<"Beginning Event Loop\n";

	//TFile * file = new TFile(in_name);
	TChain * file = new TChain("ePi");

	if( data_sim==0 ){
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim_10637_rich.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim*");

		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_*.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.4/run_skim_*.root");
		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.6/run_skim_*.root");

		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root");

		//file->Add("../trees/final_skims/10.4/tight_skim.root");
		//file->Add("../trees/final_skims/10.6/tight_skim.root");
	}
	else{
		file->Add("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim_10637.root");
	}
	
	
		//TTreeReader reader("ePi", file);
	TTreeReader reader( file);

	TTreeReaderValue<electron> e(reader, "e");

	TTreeReaderArray<pion> pi(reader, "pi");
	//TTreeReaderArray<pion> k(reader, "ka");

	TTreeReaderArray<int> pidRICH_pi(reader, "pidRICH");
	//TTreeReaderArray<int> pidRICH_k(reader, "pidRICH_ka");
	
	TTreeReaderArray<double> probRICH_pi(reader, "probRICH");
	//TTreeReaderArray<double> probRICH_k(reader, "probRICH_ka");

	TTreeReaderArray<double> nRICH_pi(reader, "nRICH");
	//TTreeReaderArray<double> nRICH_k(reader, "nRICH_ka");
	//TTreeReaderArray<bool> isGoodPion_vec(reader, "isGoodPion");
	//TTreeReaderArray<bool> isGoodPion_3d_vec(reader, "isGoodPion_3d");
	//TTreeReaderArray<bool> isGoodPion_vec_no_acc(reader, "isGoodPion_no_acc");

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);

	anal.setAnalyzerLevel(0);//runType);
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));
	int event_count = 0;
	while (reader.Next()) {
		if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
		event_count++;
		//if( event_count == 20000000 ) break;

		double p_e = e->get3Momentum().Mag();
		double theta_e = e->get3Momentum().Theta()*rad_to_deg;
		double phi_e = e->get3Momentum().Phi()* rad_to_deg;

		//if (anal.applyAcceptanceMap(p_e, phi_e, theta_e, 0) < 0) continue;
		//if( !anal.applyElectronKinematicCuts(*e) ) continue;
		//int chargeIdx = 0;
		//if( type == 0 &&  (int) ( pi.end() - pi.begin() ) == 2
		//		&& Mx_2pi_cut > 0 && 
		//		checkDiffractive(pi[0], pi[1], *e, Mx_2pi_cut)) continue;

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			//if( (matchType==2 && !isGoodPion_vec[i])||(matchType==0 && !isGoodPion_vec_no_acc[i]) ){continue;}
			if( data_sim==0 && pi[i].getBeta_rich() < .0001){continue;}
			if( data_sim==0 && nRICH_pi[i] <= nPhotons-1 ) continue;
			if( data_sim==0 && probRICH_pi[i] < prob ) continue;
			//if( !anal.applyPionKinematicCuts((pi)[i]) ) continue;

			int chargeIdx = (int)(pi[i].getCharge() < 1);
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()* rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()* rad_to_deg;
			if( p_pi < 1.25 || p_pi > 5 ) continue;
            //if (anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, chargeIdx+1) < 0) continue;

			double pid_true = (data_sim == 0)?abs(pidRICH_pi[i]):abs(pi[i].getPID());

			if(pid_true == 211 && abs(pi[i].getPID_eb()) != 2212 ){
				hP_true[chargeIdx]->Fill(p_pi);
				hTheta_p_true[chargeIdx]->Fill(p_pi, theta_pi);
				hNRICH_p[chargeIdx]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx]->Fill(p_pi, probRICH_pi[i]);
				
				if(  abs(pi[i].getPID_eb()) == 211 ){
					hP_eb[chargeIdx]->Fill(p_pi);
					hTheta_p_eb[chargeIdx]->Fill(p_pi, theta_pi);
					hNRICH_p_eb[chargeIdx]->Fill(p_pi, nRICH_pi[i]);
					hProbRICH_p_eb[chargeIdx]->Fill(p_pi, probRICH_pi[i]);
				}
			}

			if( pid_true == 321 && abs(pi[i].getPID_eb()) != 2212 ){
				hP_true[chargeIdx+2]->Fill(p_pi);
				hTheta_p_true[chargeIdx+2]->Fill(p_pi, theta_pi);
				
				hNRICH_p[chargeIdx+2]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx+2]->Fill(p_pi, probRICH_pi[i]);
				if( abs(pi[i].getPID_eb()) == 321 ){
					hP_eb[chargeIdx+2]->Fill(p_pi);
					hTheta_p_eb[chargeIdx+2]->Fill(p_pi, theta_pi);
					hNRICH_p_eb[chargeIdx+2]->Fill(p_pi, nRICH_pi[i]);
					hProbRICH_p_eb[chargeIdx+2]->Fill(p_pi, probRICH_pi[i]);
				}
			}

			if( pid_true == 2212 ){
				hP_true[chargeIdx+4]->Fill(p_pi);
				hTheta_p_true[chargeIdx+4]->Fill(p_pi, theta_pi);
				
				hNRICH_p[chargeIdx+4]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx+4]->Fill(p_pi, probRICH_pi[i]);
				if( abs(pi[i].getPID_eb()) == 2212 ){
					hP_eb[chargeIdx+4]->Fill(p_pi);
					hTheta_p_eb[chargeIdx+4]->Fill(p_pi, theta_pi);
				}
			}

		}

	}


	outFile->cd();
	
	
	for( int i = 0; i < 6; i++ ){
		hP_true[i]->Write();
		hP_eb[i]->Write();
		hTheta_p_eb[i]->Write();
		hTheta_p_true[i]->Write();
		//hP_opp[i]->Write();
		hNRICH_p[i]->Write();
		hProbRICH_p[i]->Write();
		hNRICH_p_eb[i]->Write();
		hProbRICH_p_eb[i]->Write();
	}

	// Ratios: pip/pim (indices 0,1) and kp/km (indices 2,3)
	TString ratio_species[6] = {"pip","pim", "kp", "km", "prtp", "prtm"};
	for( int i = 0; i < 6; i++ ){
		TH1F * hRatio_p = (TH1F*) hP_eb[i]->Clone("hRatio_p_" + ratio_species[i]);
		hRatio_p->Divide(hP_true[i]);
		hRatio_p->SetTitle("p^{+}/p^{-} ratio (" + ratio_species[i] + ");p (GeV/c);N^{+}/N^{-}");
		hRatio_p->Write();

		TH2F * hRatio_theta_p = (TH2F*) hTheta_p_eb[i]->Clone("hRatio_theta_p_" + ratio_species[i]);
		hRatio_theta_p->Divide(hTheta_p_true[i]);
		hRatio_theta_p->SetTitle("#theta vs p ratio (" + ratio_species[i] + ");p (GeV/c);#theta (deg)");
		hRatio_theta_p->Write();
	}
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}


bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut ){
	TLorentzVector missing = p_rest + e.getQ() - pi_1.get4Momentum() - pi_2.get4Momentum();
	return missing.Mag2() < Mx_cut;
}

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "analyzer.h"
#include "DCfid_SIDIS.h"
#define HIST_PATH _HIST

using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;

bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut );


int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File]\n";
		cerr << "[Type (0 - pi2k, 1 - k2pi)] [Theta Cut (optional)] [Mx_2pi^2 cut (optional, default 1.25)]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";


	int p_bins = 20;

    TString out_name = argv[1];
		
	int nPhotons = atoi(argv[2]);
	double prob = atof(argv[3]);
	int data_sim = atoi(argv[4]);

    TFile * outFile = new TFile(out_name, "RECREATE");//(TString) HIST_PATH + "/" + out_name + ".root", "RECREATE");
	
	// Declare histograms

	cout<<"Creating Histograms\n";

	// 4 momentum histograms — one per species
	// Index: 0=pip misid, 1=pim misid, 2=kp, 3=km
	TString species[6] = {"pip", "pim", "kp", "km", "prot", "aprot"};
	TH1F * hP_eb[6];
	TH1F * hP_true[6];
	TH1F * hP_opp[6];


	TH2F * hTheta_p_eb[6];
	TH2F * hTheta_p_true[6];
	TH2F * hNRICH_p[6];
	TH2F * hProbRICH_p[6];
	TH2F * hNRICH_p_eb[6];
	TH2F * hProbRICH_p_eb[6];

	for( int i = 0; i < 6; i++ ){
		hP_eb[i] = new TH1F(
			"hP_eb_" + species[i],
			"p (" + species[i] + ");p (GeV/c);Counts",
			p_bins, 1.25, 5.0 );
		hP_true[i] = new TH1F(
			"hP_true_" + species[i],
			"p (" + species[i] + ");p (GeV/c);Counts",
			p_bins, 1.25, 5.0 );
		hTheta_p_eb[i] = new TH2F(
			"hTheta_p_eb_" + species[i],
			"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
			50, 1.25, 5.0, 30, 5.0, 35.0 );
		hTheta_p_true[i] = new TH2F(
			"hTheta_p_true_" + species[i],
			"#theta vs p (" + species[i] + ");p (GeV/c);#theta (deg)",
			50, 1.25, 5.0, 30, 5.0, 35.0 );
		hP_opp[i] = new TH1F(
				"hP_opp_" + species[i],
				"p (" + species[i] + ");p (GeV/c);Counts",
				p_bins, 1.25, 5.0 );
		hNRICH_p[i] = new TH2F(
				"hNRICH_p_" + species[i],
				"n_{RICH} vs p (" + species[i] + ");p (GeV/c);n_{RICH}",
				50, 1.25, 5.0, 50, 0, 50 );
		hProbRICH_p[i] = new TH2F(
				"hProbRICH_p_" + species[i],
				"prob_{RICH} vs p (" + species[i] + ");p (GeV/c);prob_{RICH}",
				50, 1.25, 5.0, 50, 0, 1.0 );
		hNRICH_p_eb[i] = new TH2F(
				"hNRICH_p_eb_" + species[i],
				"n_{RICH} vs p, RICH-matched (" + species[i] + ");p (GeV/c);n_{RICH}",
				50, 1.25, 5.0, 50, 0, 50 );
		hProbRICH_p_eb[i] = new TH2F(
				"hProbRICH_p_eb_" + species[i],
				"prob_{RICH} vs p, RICH-matched (" + species[i] + ");p (GeV/c);prob_{RICH}",
				50, 1.25, 5.0, 50, 0, 1.0 );
	}

	cout<<"Beginning Event Loop\n";

	//TFile * file = new TFile(in_name);
	TChain * file = new TChain("ePi");

	if( data_sim==0 ){
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim_10637_rich.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim*");

		file->Add("/volatile/clas12/users/jphelan/SIDIS/data/detector_skims/10.2/run_skim_*.root");

		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.2/final_skim.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.4/final_skim.root");
		//file->Add("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/10.6/final_skim.root");

		//file->Add("../trees/final_skims/10.4/tight_skim.root");
		//file->Add("../trees/final_skims/10.6/tight_skim.root");
	}
	else{
		file->Add("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/detector_skim_10637_rich.root");
	}
	
	
		//TTreeReader reader("ePi", file);
	TTreeReader reader( file);

	TTreeReaderValue<electron> e(reader, "e");

	TTreeReaderArray<pion> pi(reader, "pi");
	//TTreeReaderArray<pion> k(reader, "ka");

	TTreeReaderArray<int> pidRICH_pi(reader, "pidRICH");
	//TTreeReaderArray<int> pidRICH_k(reader, "pidRICH_ka");
	
	TTreeReaderArray<double> probRICH_pi(reader, "probRICH");
	//TTreeReaderArray<double> probRICH_k(reader, "probRICH_ka");

	TTreeReaderArray<double> nRICH_pi(reader, "nRICH");
	//TTreeReaderArray<double> nRICH_k(reader, "nRICH_ka");
	//TTreeReaderArray<bool> isGoodPion_vec(reader, "isGoodPion");
	//TTreeReaderArray<bool> isGoodPion_3d_vec(reader, "isGoodPion_3d");
	//TTreeReaderArray<bool> isGoodPion_vec_no_acc(reader, "isGoodPion_no_acc");

	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);

	anal.setAnalyzerLevel(0);//runType);
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));
	int event_count = 0;
	while (reader.Next()) {
		if(event_count%100000 == 0){cout<<"Events Analyzed: "<<event_count<<std::endl;}
		event_count++;
		if( event_count == 20000000 ) break;

		double p_e = e->get3Momentum().Mag();
		double theta_e = e->get3Momentum().Theta()*rad_to_deg;
		double phi_e = e->get3Momentum().Phi()* rad_to_deg;

		//if (anal.applyAcceptanceMap(p_e, phi_e, theta_e, 0) < 0) continue;
		//if( !anal.applyElectronKinematicCuts(*e) ) continue;
		//int chargeIdx = 0;
		//if( type == 0 &&  (int) ( pi.end() - pi.begin() ) == 2
		//		&& Mx_2pi_cut > 0 && 
		//		checkDiffractive(pi[0], pi[1], *e, Mx_2pi_cut)) continue;

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			//if( (matchType==2 && !isGoodPion_vec[i])||(matchType==0 && !isGoodPion_vec_no_acc[i]) ){continue;}
			if( data_sim==0 && pi[i].getBeta_rich() < .0001){continue;}
			if( data_sim==0 && nRICH_pi[i] <= nPhotons-1 ) continue;
			if( data_sim==0 && probRICH_pi[i] < prob ) continue;
			//if( !anal.applyPionKinematicCuts((pi)[i]) ) continue;

			int chargeIdx = (int)(pi[i].getCharge() < 1);
			double p_pi = pi[i].get3Momentum().Mag();
			double theta_pi = pi[i].get3Momentum().Theta()* rad_to_deg;
			double phi_pi = pi[i].get3Momentum().Phi()* rad_to_deg;
			if( p_pi < 1.25 || p_pi > 5 ) continue;
            //if (anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, chargeIdx+1) < 0) continue;

			double pid_true = (data_sim == 0)?abs(pidRICH_pi[i]):abs(pi[i].getPID());

			if(pid_true == 211 ){
				hP_true[chargeIdx]->Fill(p_pi);
				hTheta_p_true[chargeIdx]->Fill(p_pi, theta_pi);
				hNRICH_p[chargeIdx]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx]->Fill(p_pi, probRICH_pi[i]);
				
				if(  abs(pi[i].getPID_eb()) == 211 ){
					hP_eb[chargeIdx]->Fill(p_pi);
					hTheta_p_eb[chargeIdx]->Fill(p_pi, theta_pi);
					hNRICH_p_eb[chargeIdx]->Fill(p_pi, nRICH_pi[i]);
					hProbRICH_p_eb[chargeIdx]->Fill(p_pi, probRICH_pi[i]);
				}
			}

			if( pid_true == 321 ){
				hP_true[chargeIdx+2]->Fill(p_pi);
				hTheta_p_true[chargeIdx+2]->Fill(p_pi, theta_pi);
				
				hNRICH_p[chargeIdx+2]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx+2]->Fill(p_pi, probRICH_pi[i]);
				if( abs(pi[i].getPID_eb()) == 321 ){
					hP_eb[chargeIdx+2]->Fill(p_pi);
					hTheta_p_eb[chargeIdx+2]->Fill(p_pi, theta_pi);
					hNRICH_p_eb[chargeIdx+2]->Fill(p_pi, nRICH_pi[i]);
					hProbRICH_p_eb[chargeIdx+2]->Fill(p_pi, probRICH_pi[i]);
				}
			}

			if( pid_true == 2212 ){
				hP_true[chargeIdx+4]->Fill(p_pi);
				hTheta_p_true[chargeIdx+4]->Fill(p_pi, theta_pi);
				
				hNRICH_p[chargeIdx+4]->Fill(p_pi, nRICH_pi[i]);
				hProbRICH_p[chargeIdx+4]->Fill(p_pi, probRICH_pi[i]);
				if( abs(pi[i].getPID_eb()) == 2212 ){
					hP_eb[chargeIdx+4]->Fill(p_pi);
					hTheta_p_eb[chargeIdx+4]->Fill(p_pi, theta_pi);
				}
			}

		}

	}


	outFile->cd();
	
	
	for( int i = 0; i < 6; i++ ){
		hP_true[i]->Write();
		hP_eb[i]->Write();
		hTheta_p_eb[i]->Write();
		hTheta_p_true[i]->Write();
		//hP_opp[i]->Write();
		hNRICH_p[i]->Write();
		hProbRICH_p[i]->Write();
		hNRICH_p_eb[i]->Write();
		hProbRICH_p_eb[i]->Write();
	}

	// Ratios: pip/pim (indices 0,1) and kp/km (indices 2,3)
	TString ratio_species[6] = {"pip","pim", "kp", "km", "prtp", "prtm"};
	for( int i = 0; i < 6; i++ ){
		TH1F * hRatio_p = (TH1F*) hP_eb[i]->Clone("hRatio_p_" + ratio_species[i]);
		hRatio_p->Divide(hP_true[i]);
		hRatio_p->SetTitle("p^{+}/p^{-} ratio (" + ratio_species[i] + ");p (GeV/c);N^{+}/N^{-}");
		hRatio_p->Write();

		TH2F * hRatio_theta_p = (TH2F*) hTheta_p_eb[i]->Clone("hRatio_theta_p_" + ratio_species[i]);
		hRatio_theta_p->Divide(hTheta_p_true[i]);
		hRatio_theta_p->SetTitle("#theta vs p ratio (" + ratio_species[i] + ");p (GeV/c);#theta (deg)");
		hRatio_theta_p->Write();
	}
	cout<<"FINISHED WRITING\n";
	outFile->Close();
}


bool checkDiffractive( pion pi_1, pion pi_2, electron e, double Mx_cut ){
	TLorentzVector missing = p_rest + e.getQ() - pi_1.get4Momentum() - pi_2.get4Momentum();
	return missing.Mag2() < Mx_cut;
}