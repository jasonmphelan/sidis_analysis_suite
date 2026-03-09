#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
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
#include "analyzer.h"
#include "constants.h"
#include "cut_values.h"

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;
using std::isnan;
using namespace cutVals;
using namespace constants;

double getVarVal(TString var,  electron e, pion pi );

int main( int argc, char** argv){

	if( argc <3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File]\n";
		return -1;
	}
	cerr << "Files used: " << argv[1] << " " << argv[2] <<"\n";

	TString in_name = argv[1];
    TString out_name = argv[2];
	double eBeam = atof(argv[3]);
	double Mx_low = atof(argv[4]);
	double Mx_high = atof(argv[5]);
	double norm_bound = atof(argv[6]);
	TString var_name = argv[7];
	int bins_var = atoi(argv[8]);
	double var_min = atof(argv[9]);
	double var_max = atof(argv[10]);


	analyzer anal( 0, -1 );
	anal.setAnalyzerLevel(0);//runType);
	anal.loadMatchingFunctions("matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");//%.1f.root", energy));


	int nBinsXb = 14;
	int nBinsQ2 = 12;
	int nBinsZ = 14;

	double widthXb = (.66 - .1)/( (double) nBinsXb );
	double widthQ2 = (8. - 2.)/( (double) nBinsQ2 );
	double widthZ = (1. - .3)/( (double) nBinsZ );

        //TFile * file_rec = new TFile(in_name);
	TFile * file_out = new TFile(out_name + in_name + ".root", "RECREATE");
	TVector3 * bounds = new TVector3(Mx_low, Mx_high, norm_bound);

	TH1F* hMx_2pi[2][bins_var][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho[2][bins_var][nBinsXb][nBinsQ2][nBinsZ];
	TH1F* hM_rho_back[2][bins_var][nBinsXb][nBinsQ2][nBinsZ];


	TH3F * hNorms[bins_var][2];
	for( int var = 0; var < bins_var; var ++ ){
		hNorms[var][0] = new TH3F( Form("hNorm_pip_%i", var), "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);
		hNorms[var][1] = new TH3F( Form("hNorm_pim_%i", var), "", nBinsXb, .1, .66, nBinsQ2, 2, 8, nBinsZ, .3, 1);

		for( int x = 0; x < nBinsXb; x++ ){
			for( int q = 0; q < nBinsQ2; q++ ){
				for( int z = 0; z <  nBinsZ; z++ ){
					TString binLims = Form("%.2f<xB<%.2f, %.1f<Q2<%.1f, %.2f<Z<%.2f", .1 + x*widthXb, .1 + (x+1)*widthXb, 2 + q*widthQ2, 2 + (q+1)*widthQ2, .3 + z*widthZ, .3 + (z + 1)*widthZ);
					hMx_2pi[0][var][x][q][z] = new TH1F( Form("hMx_pip_%i_%i_%i_%i", var, x, q, z), "pip " + binLims, 75, 0, 2.5);
					hMx_2pi[1][var][x][q][z] = new TH1F( Form("hMx_pim_%i_%i_%i_%i", var, x, q, z), "pim " + binLims, 75, 0, 2.5);
					
					hM_rho[0][var][x][q][z] = new TH1F( Form("hM_rho_pip_%i_%i_%i_%i", var, x, q, z), "pip " + binLims, 75, 0, 2.5);
					hM_rho[1][var][x][q][z] = new TH1F( Form("hM_rho_pim_%i_%i_%i_%i", var, x, q, z), "pim " + binLims, 75, 0, 2.5);
					hM_rho_back[0][var][x][q][z] = new TH1F( Form("hM_rho_back_pip_%i_%i_%i_%i", var, x, q, z), "pip " + binLims, 75, 0, 2.5);
					hM_rho_back[1][var][x][q][z] = new TH1F( Form("hM_rho_back_pim_%i_%i_%i_%i", var, x, q, z), "pim " + binLims, 75, 0, 2.5);
				}
			}
		}
	}
	TChain * file_rec = new TChain("ePi");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.2" + in_name + ".root");
	if(eBeam == 10.2 || eBeam == 0 )file_rec->Add("../trees/final_skims/rho_skims/rotated_10.2" + in_name + ".root");
	if(eBeam == 10.4 || eBeam == 0 )file_rec->Add("../trees/final_skims/rho_skims/rotated_10.4" + in_name + ".root");
	if(eBeam == 10.6 || eBeam == 0 )file_rec->Add("../trees/final_skims/rho_skims/rotated_10.6" + in_name + ".root");

	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.4" + in_name + ".root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rotated_10.6" + in_name + ".root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.2.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.4.root");
	//file_rec->Add("/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/rho_skim_10.6.root");


	//Load input tree
    TTreeReader reader_rec(file_rec);

	TTreeReaderValue<electron> e(reader_rec, "e");
	TTreeReaderValue<double> Mx_2pi(reader_rec, "Mx_2pi");
	TTreeReaderValue<double> M_rho(reader_rec, "M_rho");
	TTreeReaderArray<pion> pi(reader_rec, "pi");
	TTreeReaderArray<bool> isGoodPion(reader_rec, "isGoodPion"); //should be no acceptance
	TTreeReaderArray<double> rhoWeight(reader_rec, "rhoWeight");
	TTreeReaderArray<double> rhoWeight_sym(reader_rec, "rhoWeight_sym");

	//Fill histograms
	while (reader_rec.Next()) {
                int event_count = reader_rec.GetCurrentEntry();

		if(event_count%100000 == 0){
			cout<<"Events Analyzed: "<<event_count<<std::endl;
		}
	
		for( int i = 0; i < (int)(pi.end() - pi.begin());i++ ){
			if( !isGoodPion[i] ){continue;}
			//if( !isGoodPion[0] || !isGoodPion[1]){continue;}
			int chargeIdx = (int)( pi[i].getCharge() < 0 );
			int this_bin_Q2 = (int)( ( (e->getQ2() - Q2_min)/(Q2_max-Q2_min) )*nBinsQ2);
			int this_bin_xB = (int)( ( (e->getXb() - xB_min)/(xB_max-xB_min) )*nBinsXb);
			int this_bin_var = (int)( ( (getVarVal(var_name, *e, pi[i]) - var_min)/(var_max - var_min) )*bins_var);
			int this_bin_Z = (int)( ( (pi[i].getZ() - .3)/(1.-.3) )*nBinsZ);
			
			bool matching = anal.applyAcceptanceMatching(pi[i], 2);
				//matching = isGoodPion[i]; }
			
			if( !matching ){ continue; }


			if( rhoWeight[i] <= 0 || rhoWeight[i] > 100 ){ continue; }
		
			hMx_2pi[chargeIdx][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*Mx_2pi) , 0.5*rhoWeight[i] );
			if( *Mx_2pi > bounds->X() && *Mx_2pi < bounds->Y() ) hM_rho_back[chargeIdx][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) , 0.5*rhoWeight[i] );
			if( *Mx_2pi < bounds->X() ) hM_rho[chargeIdx][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) , 0.5*rhoWeight[i] );

			if( rhoWeight_sym[i] <= 0 || rhoWeight_sym[i] > 100 ){ continue; }
		
			hMx_2pi[(int)(!chargeIdx)][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*Mx_2pi) , 0.5*rhoWeight_sym[i] );
			if( *Mx_2pi > bounds->X() && *Mx_2pi < bounds->Y() ) hM_rho_back[(int)(!chargeIdx)][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) , 0.5*rhoWeight_sym[i] );
			if( *Mx_2pi < bounds->X() ) hM_rho[(int)(!chargeIdx)][this_bin_var][this_bin_xB][this_bin_Q2][this_bin_Z]->Fill( (*M_rho) , 0.5*rhoWeight_sym[i] );

		}
	}

	
	TFile * out_fits = new TFile(out_name + "_M_rho.root", "RECREATE");
	out_fits->cd();

	for( int var = 0; var < bins_var; var++ ){
		for( int x = 0; x < nBinsXb; x++ ){
			for( int q = 0; q < nBinsQ2; q++ ){
				for( int z = 0; z <  nBinsZ; z++ ){				
				
					double bin_max= hM_rho[0][var][x][q][z]->FindBin( bounds->Z() );
					double back_num_pip = 0;
					double back_num_pim = 0;
					double rho_num_pip = 0;
					double rho_num_pim = 0;

					for( int bin = 1; bin <= bin_max; bin++ ){
						back_num_pip+=hM_rho_back[0][var][x][q][z]->GetBinContent(bin);
						back_num_pim+=hM_rho_back[1][var][x][q][z]->GetBinContent(bin);
						rho_num_pip+=hM_rho[0][var][x][q][z]->GetBinContent(bin);
						rho_num_pim+=hM_rho[1][var][x][q][z]->GetBinContent(bin);
					}


					double norm_pip = rho_num_pip/back_num_pip;//fit_dis_pip->Integral( .75, 1.15 )/fit_dis_pip->Integral(1.15, 1.45);
					double norm_pim = rho_num_pim/back_num_pim;//fit_dis_pim->Integral( .75, 1.15 )/fit_dis_pim->Integral(1.15, 1.45);
					if( norm_pip < 0 || isnan(norm_pip) || !isfinite(norm_pip) ) norm_pip = 0;
					if( norm_pim < 0 || isnan(norm_pim) || !isfinite(norm_pim) ) norm_pim = 0;
					
					if( norm_pip > 0 && !isnan(norm_pip) && isfinite(norm_pip) )cout<<"Norm pip bin "<<z*0.05+.3<<" : "<<norm_pip<<std::endl;

					hNorms[var][0]->SetBinContent( x+1, q+1, z+1, norm_pip );
					hNorms[var][1]->SetBinContent( x+1, q+1, z+1, norm_pim );
					
					hMx_2pi[0][var][x][q][z]->Write();
					hMx_2pi[1][var][x][q][z]->Write();
				}
			}
		}
	}	
	out_fits->Close();
	
	file_out->cd();
	bounds->Write("bounds");
	for( int var = 0; var < bins_var; var++ ){
		hNorms[var][0]->Write();
		hNorms[var][1]->Write();
	}
	file_out->Close();

}

double getVarVal(TString var,  electron e, pion pi ){
	if( var == "p_e" ) return e.get3Momentum().Mag();
	if( var == "theta_e" ) return e.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_e" ) return e.get3Momentum().Phi()*rad_to_deg;
	if( var == "W2" ) return e.getW2();
	if( var == "Q2" ) return e.getQ2();
	if( var == "xB" ) return e.getXb();
	if( var == "y" ) return e.getY();
	if( var == "sector" ) return e.getDC_sector();
	if( var == "sector_e" ) { double phi_deg = e.get3Momentum().Phi()*rad_to_deg + 20.; if (phi_deg < 0.) phi_deg += 360.; return (int)(phi_deg / 60.) + 1; }


	if( var == "p_pi" ) return pi.get3Momentum().Mag();
	if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
	if( var == "phi_pi" ) return pi.get3Momentum().Phi()*rad_to_deg;
	if( var == "phi_q" ) return pi.getPi_q().Phi()*rad_to_deg;
	if( var == "Z" || var == "z" ) return pi.getZ();
	if( var == "Mx" || var == "M_x" ) return pi.getMx();
	if( var == "pT" || var == "Pt" ) return  pi.getPi_q().Pt();

	return 0;
}
