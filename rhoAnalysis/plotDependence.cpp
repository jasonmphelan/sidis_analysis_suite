
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
#include "correctionTools.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST


using std::cerr;
using std::isfinite;
using std::cout;
using std::ofstream;

using namespace cutVals;
using namespace constants;


void zeroSuppress( TH1F * h);
double getVarVal( electron e, pion pi ); 

int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [input r file] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " <<(TString) HIST_PATH +"/" + argv[2] <<"\n";

	TString in_name = argv[1];
       	TString r_name = argv[2];
	TString out_name = argv[3];
	int matchType = 2;
	int applyCorr = 4;

	cout<<"Create and Load files\n";

	TFile * outFile = new TFile((TString) HIST_PATH + "/" + out_name, "RECREATE");

	//TFile * outFile = new TFile( outName, "RECREATE");
	TFile * inFile = new TFile( in_name );
	TFile * rFile = new TFile( r_name );

	double sumWeights[4][bins_Q2+1][bins_xB+1][bins_Z + 1][2] = {0};
	double sumWeightsErr[4][bins_Q2+1][bins_xB+1][bins_Z + 1][2] = {0};

	TH1F * hZ[bins_Q2+1][bins_xB+1][2];
	TH1F * hZ_r[bins_Q2+1][bins_xB+1][2];
	TH1F * hZ_r_bac[bins_Q2+1][bins_xB+1][2];
	
	TH1F * hQ2[bins_xB+1][2];
	TH1F * hQ2_r[bins_xB+1][2];
	TH1F * hQ2_r_bac[bins_xB+1][2];

	TString charge_str[2] = {"_Pip", "_Pim"};

	cout<<"Make histograms\n";

	for( int i = 0; i < bins_Q2; i++ ){
		for( int j = 0; j < bins_xB; j++ ){
			for( int k = 0; k < 2; k++ ){
			
				TString title = Form("%.1f < Q^{2} < %.1f [GeV^{2}], %.2f < x_{B} < %.2f", 
							Q2_min + (double)i*(0.5), Q2_min + (double)(i+1.)*(0.5), 
							xB_min + (double)j*(0.04), xB_min + (double)(j+1.)*(0.04)); 

				hZ[i][j][k] = new TH1F("hZ" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , title , bins_Z, .3, 1);
				hZ_r[i][j][k] = new TH1F("hZ_r" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , title , bins_Z, .3, 1);
				hZ_r_bac[i][j][k] = new TH1F("hZ_r_bac" + charge_str[k] + Form("_%i_%i",  i+1, j+1) , title , bins_Z, .3, 1);
						
			}
		}
	}
	for( int j = 0; j < bins_xB; j++ ){
		for( int k = 0; k < 2; k++ ){
		
			TString title = Form("%.2f < x_{B} < %.2f",  
						xB_min + (double)j*(0.04), xB_min + (double)(j+1.)*(0.04)); 

			hQ2[j][k] = new TH1F("hQ2" + charge_str[k] + Form("_%i", j+1) , title , bins_Z, .3, 1);
			hQ2_r[j][k] = new TH1F("hQ2_r" + charge_str[k] + Form("_%i", j+1) , title , bins_Z, .3, 1);
			hQ2_r_bac[j][k] = new TH1F("hQ2_r_bac" + charge_str[k] + Form("_%i",  j+1) , title , bins_Z, .3, 1);
						
		}
	}
	
	correctionTools corrector(1);
	corrector.loadFits();	

	////////////////////////////
	///////// Pions ////////////
	////////////////////////////

	cout<<"Begin analysis\n";
	double beam_energy = 10.2; //current energy of file

	TTreeReader reader_rec("ePi", inFile);

	TTreeReaderValue<TLorentzVector> beam( reader_rec, "beam" );
	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion");

	int event_total = reader_rec.GetEntries();
	double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	while (reader_rec.Next()) {
		int event_count = reader_rec.GetCurrentEntry();
		
		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
		}

		if( beam->E() != beam_energy ){
			corrector.loadNewEnergy( beam->E() );
			beam_energy = beam->E();
		}

		for( int i = 0; i < (int) ( pi.end() - pi.begin() ); i++ ){
			
			int chargeIdx = (int)( pi[i].getCharge() < 1 );
			double p_pi = pi[i].get3Momentum().Mag();
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*bins_Z);

			bool matching = true;
			if( matchType == 2 ){ matching = !isGoodPion[i]; }
			//else if( matchType == 3 ){ matching = !isGoodPion3d[i]; }
			//else{ matching = false; }

			if( matching ){ continue; }
			corrector.setKinematics( e->getXb(), e->getQ2(), pi[i].getZ(), p_pi );

			double weight = 1;
			
			//double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
			//double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
			double mc_weight = corrector.getCorrectionFactor( 1, chargeIdx );
			double k_weight = corrector.getCorrectionFactor(2, chargeIdx);
			//double bin_err = corrector.getCorrectionError(0, chargeIdx);
			//double acc_err = corrector.getCorrectionError(1, chargeIdx);
			double mc_err = corrector.getCorrectionFactor( 1, chargeIdx );
			double k_err = corrector.getCorrectionError(2, chargeIdx);
			
			if( applyCorr > 0 ){
				//if( !isfinite(acc_weight) || acc_err/acc_weight > .2){acc_weight = 0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
				//if( !isfinite(bin_weight) || bin_err/bin_weight > .2){bin_weight = 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
				
				//weight *= acc_weight*bin_weight;
				weight *= mc_weight;
			}
			if( applyCorr > 1 ){
				weight *= k_weight;
			}	

			double eventWeightErr = 0;
			if( applyCorr > 1 ){
				//eventWeightErr +=  weight*sqrt( pow(bin_err/bin_weight, 2) + pow(acc_err/acc_weight, 2) + pow(k_err/k_weight, 2) );
				eventWeightErr +=  weight*sqrt( pow(mc_err/mc_weight, 2) + pow(k_err/k_weight, 2) );
			}
			if( applyCorr == 1 ){
				//eventWeightErr += weight*sqrt( pow(bin_err/bin_weight, 2) + pow(acc_err/acc_weight, 2) );
				eventWeightErr += mc_err;
			}

			//weight *= corrector.getCorrectionFactor(2, chargeIdx);	
			
			sumWeights[0][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += weight;	
			sumWeightsErr[0][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += eventWeightErr;	

			hZ[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( pi[i].getZ(), weight );
			hQ2[this_bin_xB][chargeIdx]->Fill( e->getQ2(), weight );
			//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
			

		}
	}
	
	
	///////////////////////////////////////
	////////////// Rhos /////////////////
	///////////////////////////////////////
	
	TTreeReader reader_r("ePi", rFile);
	TTreeReaderValue<TLorentzVector> beam_r( reader_r, "beam" );

	TTreeReaderValue<electron> e_r(reader_r, "e");
	TTreeReaderArray<pion> r(reader_r, "pi");
	TTreeReaderValue<double> rhoWeight( reader_r, "rhoWeight");
	TTreeReaderValue<double> Mx_2pi( reader_r, "Mx_2pi");
	TTreeReaderValue<double> rhoError( reader_r, "rhoErr");
	TTreeReaderArray<bool> isGoodRho(reader_r, "isGoodPion");

	event_total = reader_r.GetEntries();
	//double events_in_bin[2][bins_Q2][bins_xB][bins_Z][bins_p] = {0};

	if( applyCorr>-1 ){
		while (reader_r.Next()) {
			int event_count = reader_r.GetCurrentEntry();

			if(event_count%10000 == 0){
				cout<<"Events Analyzed: "<<event_count<< " / "<<event_total<<std::endl;
			}
			if( beam_r->E() != beam_energy ){
				std::cout<<"Setting correction energy\n";
				corrector.loadNewEnergy( beam_r->E() );
				beam_energy = beam_r->E();
			}

			for( int i = 0; i < (int) ( r.end() - r.begin() ); i++ ){
				
				int chargeIdx = (int)( r[i].getCharge() < 1 );
				double p_pi = r[i].get3Momentum().Mag();
				int this_bin_Q2 = (int)( ( (e_r->getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
				int this_bin_xB = (int)( ( (e_r->getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);
				int this_bin_Z = (int)( ( (r[i].getZ() - .3)/(1.-.3) )*bins_Z);

				bool matching = true;
				if( matchType == 2 ){ matching = !isGoodRho[i]; }

				if( matching ){ continue; }

				corrector.setKinematics( e_r->getXb(), e_r->getQ2(), r[i].getZ(), p_pi );
				double weight = 1;
				
				//double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
				//double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
				//double bin_err = corrector.getCorrectionError(0, chargeIdx);
				//double acc_err = corrector.getCorrectionError(1, chargeIdx);
				
				//double bin_weight = corrector.getCorrectionFactor(0, chargeIdx);
				//double acc_weight = corrector.getCorrectionFactor(1, chargeIdx);
				double mc_weight = corrector.getCorrectionFactor( 1, chargeIdx );
				//double k_weight = corrector.getCorrectionFactor(3, chargeIdx);
				//double bin_err = corrector.getCorrectionError(0, chargeIdx);
				//double acc_err = corrector.getCorrectionError(1, chargeIdx);
				double mc_err = corrector.getCorrectionFactor( 1, chargeIdx );
				//double k_err = corrector.getCorrectionError(3, chargeIdx);
				
				if( applyCorr > 0 ){
					//if( !isfinite(acc_weight) || acc_err/acc_weight > .2){acc_weight = 0;}// || acc_weight < 0.2 || acc_weight > 6 ){continue;}
					//if( !isfinite(bin_weight) || bin_err/bin_weight > .2){bin_weight = 0;}// || bin_weight < 0.2 || bin_weight > 3 ){continue;}
					
					weight *= mc_weight;//acc_weight*bin_weight;
				}
				if( applyCorr > 3 ){
					if( *rhoWeight <= 1 || *rhoWeight > 10 ){continue;}
					weight *= (*rhoWeight);
				}	
					
				double eventWeightErr = 0;
				
				//eventWeightErr += weight*sqrt( pow(bin_err/bin_weight, 2) + pow(acc_err/acc_weight, 2) + pow(*rhoError / *rhoWeight, 2) );
				eventWeightErr += weight*sqrt( pow(mc_err/mc_weight, 2) + pow(*rhoError / *rhoWeight, 2) );
		
				if( isnan(eventWeightErr) ){continue;}

				if( *Mx_2pi < 1.15 ){
					hZ_r[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( r[i].getZ(), weight );
					hQ2_r[this_bin_xB][chargeIdx]->Fill( e_r->getQ2(), weight );
					sumWeights[2][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += weight*weight;	
					sumWeightsErr[2][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += eventWeightErr*eventWeightErr;	
				}
				if( *Mx_2pi > 1.15 && *Mx_2pi < 1.45 ){
					hZ_r_bac[this_bin_Q2][this_bin_xB][chargeIdx]->Fill( r[i].getZ(), weight );
					hQ2_r_bac[this_bin_xB][chargeIdx]->Fill( e_r->getQ2(), weight );
					sumWeights[3][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += weight*weight;	
					sumWeightsErr[3][this_bin_Q2][this_bin_xB][this_bin_Z][chargeIdx] += eventWeightErr*eventWeightErr;	
				}
				//events_in_bin[chargeIdx][this_bin_Q2][this_bin_xB][this_bin_Z][this_bin_p]++;
				

			}
		}
	}




	//////////////////////////////////////
	///////////// Calcs //////////////////
	//////////////////////////////////////

	cout<<"Start calculations\n";

	TH1F * helper_1 = new TH1F("helper1", "helper1", 14, .3, 1);
	TH1F * helper_2 = new TH1F("helper2", "helper2", 14, .3, 1);
		
	TFile * rho_norms = new TFile("/work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/rho_norms.root");
	TH3F * hNorms_pip = (TH3F*)rho_norms->Get("hNorm_pip");
	TH3F * hNorms_pim = (TH3F*)rho_norms->Get("hNorm_pim");
	cout<<"Found rho norms\n";
	outFile->cd();
	for( int i = 1; i <= bins_Z; i++ ){
		helper_1->SetBinContent(i, -1.);
		helper_2->SetBinContent(i, 4.);
		
		helper_1->SetBinError(i, 0.);
		helper_2->SetBinError(i, 0.);
	}

	for( int i = 1; i <= bins_Q2; i++ ){
		for( int j = 1; j <= bins_xB; j++ ){
			for( int k = 1; k <= bins_Z; k++ ){
				hZ[i-1][j-1][0]->SetBinError( k, sqrt( sumWeights[0][i-1][j-1][k-1][0] + sumWeightsErr[0][i-1][j-1][k-1][0]) );
				hZ[i-1][j-1][1]->SetBinError( k, sqrt( sumWeights[0][i-1][j-1][k-1][1] + sumWeightsErr[0][i-1][j-1][k-1][1]) );
				
				hZ_r[i-1][j-1][0]->SetBinError( k, sqrt( sumWeights[2][i-1][j-1][k-1][0] + sumWeightsErr[2][i-1][j-1][k-1][0]) );
				hZ_r[i-1][j-1][1]->SetBinError( k, sqrt( sumWeights[2][i-1][j-1][k-1][1] + sumWeightsErr[2][i-1][j-1][k-1][1]) );
				
				hZ_r_bac[i-1][j-1][0]->SetBinError( k, sqrt( sumWeights[3][i-1][j-1][k-1][0] + sumWeightsErr[3][i-1][j-1][k-1][0]) );
				hZ_r_bac[i-1][j-1][1]->SetBinError( k, sqrt( sumWeights[3][i-1][j-1][k-1][1] + sumWeightsErr[3][i-1][j-1][k-1][1]) );

			}

			cout<<"Subtracting rho\n";
			if( applyCorr == 4 ){
				for( int k = 1; k <= bins_Z; k++ ){
					double pip_cont = hZ_r_bac[i-1][j-1][0]->GetBinContent(k);
					double pip_err = hZ_r_bac[i-1][j-1][0]->GetBinError(k);
					double pip_scale = hNorms_pip->GetBinContent(j, i, k);
					double pip_scale_err = 0; //hNorms_pip->GetBinError(j, i, k);

					hZ_r_bac[i-1][j-1][0]->SetBinContent( k, pip_cont*pip_scale );
					hZ_r_bac[i-1][j-1][0]->SetBinError( k, pip_cont*pip_scale*sqrt( pow( pip_err/pip_cont, 2) 
											+ pow( pip_scale_err/pip_scale, 2) ) );
					
					double pim_cont = hZ_r_bac[i-1][j-1][1]->GetBinContent(k);
					double pim_err = hZ_r_bac[i-1][j-1][1]->GetBinError(k);
					double pim_scale = hNorms_pim->GetBinContent(j, i, k);
					double pim_scale_err = hNorms_pim->GetBinError(j, i, k);
					
					hZ_r_bac[i-1][j-1][1]->SetBinContent( k, pim_cont*pim_scale );
					hZ_r_bac[i-1][j-1][1]->SetBinError( k, pim_cont*pim_scale*sqrt( pow( pim_err/pim_cont, 2) 
											+ pow( pim_scale_err/pim_scale, 2) ) );
					

				}
			
				hZ_r[i-1][j-1][0]->Add( hZ_r_bac[i-1][j-1][0], -1 );
				hZ_r[i-1][j-1][1]->Add( hZ_r_bac[i-1][j-1][1], -1 );

				zeroSuppress(hZ_r[i-1][j-1][0]);
				zeroSuppress(hZ_r[i-1][j-1][1]);

			}

			hZ[i-1][j-1][0]->Write();
			hZ[i-1][j-1][1]->Write();
			hZ_r[i-1][j-1][0]->Write();
			hZ_r[i-1][j-1][1]->Write();
		
		}
	}		

	outFile->Close();

}


void zeroSuppress( TH1F * h){

	int nBins = h->GetNbinsX();
	for( int i = 1; i <= nBins; i++ ){
		if( h->GetBinContent(i) < 0 ){
			h->SetBinContent(i, 0);
			h->SetBinError(i, 0);
		}
		if( isnan(h->GetBinError(i)) ){
			h->SetBinError(i, 0);
		}
	}
}


double getVarVal(TString var,  electron * e, pion pi ){
	if( var == "p_e" ) return e->get3Momentum().Mag();
	if( var == "theta_e" ) return e->get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e->get3Momentum().Phi()*rad_to_deg;
	if( var == "W2" ) return e->getW2();
	if( var == "Q2" ) return e->getQ2();
	if( var == "xB" ) return e->getXb();
	if( var == "y" ) return e->getY();


	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) return pi.getPi_q().Phi()*rad_to_deg;
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	return 0;
}
