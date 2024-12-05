#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
#include "reader.h"
#include "analyzer.h"
#include "reader.h"

using namespace clas12;
using namespace constants;

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [outFile name] [Beam energy]\n";
		return -1;
	}
	
       	double Ebeam = atof(argv[2]); // [GeV]
	TString outFileName = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 && Ebeam != 0 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	anal.setAnalyzerLevel(0);
	anal.loadCutValues(-1, Ebeam);
	
	reader runReader;
	runReader.setNumFiles( 2 );
	runReader.setRunType( 0 );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
       	runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;

	//Declare hists
	TString spec[2] = {"pip", "pim"};

	TH2F * hPCAL_WV[2];
	TH2F * hEdep[2];	       
	TH2F * hSF[2][6];
	TH2F * hSF_corr[2];
	TH1F * hVz_e[2];

	TH1F * hVz_pi[2];
	TH2F * hChi2[2];
	TH2F * hFid_e[2][3][2];
	TH2F * hFid_pi[2][3][2];

	TH2F * hBeta_p_e[11];
	TH1F * hBeta_e[11];

	TH2F * hBeta_p_pi[11][2];
	TH1F * hBeta_pi[11][2];

	TH1F * hPositives = new TH1F("hPos", "hPos", 500, 0.75, 1.01);
	TH1F * hNegatives = new TH1F("hMin", "hMin", 500, 0.75, 1.01);
	
	TH2F * hPos_b_p = new TH2F("hPos_b_p", "hPos_b_p", 150, 0, 10, 500, 0.75, 1.01);
	TH2F * hNeg_b_p = new TH2F("hMin_b_p", "hMin_b_p", 150, 0, 10, 500, 0.75, 1.01);

	double e_ext[3] = {150, 300,  250};
	double pi_ext[3] = {175, 350,  275};
	TString cuts[11] = {"none", "elFid", "ePcal", "eEdep", "eSF", "eSFcorr", "eVtz", "eAll", "piFid", "piChi2", "piVtz"};
	

	for( int i = 0; i < 2; i++ ){
		hPCAL_WV[i] = new TH2F( "hPCAL_WV_"+spec[i], "", 400, 0, 200, 400, 0, 200 );
		hEdep[i] = new TH2F( "hEdep_"+spec[i], "", 400, 0, 1.4, 400, 0, 1.4 );
		hSF_corr[i] = new TH2F( "hSF_corr_"+spec[i], "", 400, 0, 2.5, 400, 0, 2.5 );
		hVz_e[i] = new TH1F( "hVz_e_"+spec[i], "", 400, -10, 10);

		hVz_pi[i] = new TH1F( "hVz_"+spec[i], "", 400, -15, 15);
		hChi2[i] = new TH2F( "hChi2_"+spec[i], "", 400, 0, 9, 400, -4.5, 4.5);
	
		for( int j = 0; j < 3; j++ ){
			hFid_e[i][j][0] = new TH2F(Form("hFid_e_bef_reg_%i_", j) + spec[i], "", 350, -e_ext[j], e_ext[j], 350, -e_ext[j], e_ext[j]);
			hFid_e[i][j][1] = new TH2F(Form("hFid_e_aft_reg_%i_", j) + spec[i], "", 350, -e_ext[j], e_ext[j], 350, -e_ext[j], e_ext[j]);
			
			hFid_pi[i][j][0] = new TH2F(Form("hFid_pi_bef_reg_%i_", j) + spec[i], "", 350, -pi_ext[j], pi_ext[j], 350, -pi_ext[j], pi_ext[j]);
			hFid_pi[i][j][1] = new TH2F(Form("hFid_pi_aft_reg_%i_", j) + spec[i], "", 350, -pi_ext[j], pi_ext[j], 350, -pi_ext[j], pi_ext[j]);
		}
		
		for( int j = 0; j < 6; j++ ){
			hSF[i][j] = new TH2F( Form("hSF_sec_%i_", j)+spec[i], "", 400, 0, 8, 400, .1, .35 );
		}
	}

	for( int i = 0; i < 11; i++ ){
		hBeta_p_e[i] = new TH2F( (TString) "hBeta_p_e_" + cuts[i], "", 400, 0, 5.5, 400, .85, 1.15);
		hBeta_p_pi[i][0] = new TH2F((TString) "hBeta_p_pip_" + cuts[i], "", 400, 0, 5.5, 400, .85, 1.15);
		hBeta_p_pi[i][1] = new TH2F((TString) "hBeta_p_pim_" + cuts[i], "", 400, 0, 5.5, 400, .85, 1.15);
	
		hBeta_e[i] = new TH1F( (TString)"hBeta_e_" + cuts[i], "", 400, .85, 1.15 );
		hBeta_pi[i][0] = new TH1F( (TString)"hBeta_pip_" + cuts[i], "", 400, .85, 1.15 );
		hBeta_pi[i][1] = new TH1F( (TString)"hBeta_pim_" + cuts[i], "", 400, .85, 1.15 );
	}

	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );

	// Electron Variables
	electron e;

	// Pion Variables
	std::vector<pion> pi ;
	
	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted
	
	// Set Output file and tree
	TFile * outputFile;
	int nFiles = 0;
	int RunType = 0;
	int inclusive = 0;
	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		//create the event reader
		clas12reader c12(files.GetFileName(i).Data());
		auto mcparts = c12.mcparts();		
		int NeventsTotal = c12.getReader().getEntries();       
    	   	int event = 0;	

    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( event%1000000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			evnum  = c12.runconfig()->getEvent();
			runnum = c12.runconfig()->getRun();
			
			///////////////////////////Initialize variables//////////////////////////////////////////////	
			electrons.clear();
			pions.clear();
			pipluses.clear();
			piminuses.clear();
			
			e.Clear();
			pi.clear();

			Ne = Npi = Npips = Npims = 0;

			/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
			//Fill hists with all negatives and positives
			std::vector<region_part_ptr> allParts = c12.getDetParticles();	
	
			for( auto particle : allParts ){
					
				if( particle->par()->getCharge() > 0 ){
					hPositives->Fill(particle->getBeta());
					hPos_b_p->Fill(particle->par()->getP(),particle->getBeta() );	
				}
				if( particle->par()->getCharge() < 0 ){
					hNegatives->Fill(particle->getBeta());
					hNeg_b_p->Fill( particle->par()->getP(),particle->getBeta());	
				}

			}

			
			// Get Particles By PID
			electrons   = c12.getByID( 11   );
			pipluses    = c12.getByID( 211  );
			piminuses   = c12.getByID(-211  );
			
			pions = pipluses;
			pions.insert( pions.end(), piminuses.begin(), piminuses.end() );

			Ne      = electrons.size();
			Npi 	= pions.size();
			
			if( Ne < 1 ){ continue; } //Keep only events with one electron...
			if( Npi == 0 && inclusive != 1 ){ continue; }	
			
			//////////////electron analysis////////////////////
			int e_idx = GetLeadingElectron(electrons, Ne);	
			e.setElectron( Ebeam, electrons[e_idx]);
			////////////////Pion analysis/////////////////
			
			pion pi_dummy;
			//Find good pions			
			for(int i = 0; i < Npi; i++){
				pi_dummy.Clear();
				pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
	
			
				int chargeIdx = (int)(pi_dummy.getCharge() < 0);
				hPCAL_WV[chargeIdx]->Fill( e.getV(), e.getW() );
				hEdep[chargeIdx]->Fill( e.getEpcal(),  e.getEecin() + e.getEecout() );
				hSF_corr[chargeIdx]->Fill( e.getEpcal()/e.get3Momentum().Mag(), (e.getEecin())/e.get3Momentum().Mag() );
				hVz_e[chargeIdx]->Fill( e.getVt().Z() );

				hVz_pi[chargeIdx]->Fill( pi_dummy.getVt().Z() - e.getVt().Z() );
				hChi2[chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getChi2() );
	
				hFid_e[chargeIdx][0][0]->Fill( e.getDC_x1(), e.getDC_y1() ); 
				hFid_pi[chargeIdx][0][0]->Fill( pi_dummy.getDC_x1(), pi_dummy.getDC_y1() );
				
				hFid_e[chargeIdx][1][0]->Fill( e.getDC_x2(), e.getDC_y2() ); 
				hFid_pi[chargeIdx][1][0]->Fill( pi_dummy.getDC_x2(), pi_dummy.getDC_y2() );
				
				hFid_e[chargeIdx][2][0]->Fill( e.getDC_x3(), e.getDC_y3() ); 
				hFid_pi[chargeIdx][2][0]->Fill( pi_dummy.getDC_x3(), pi_dummy.getDC_y3() );
				
				hSF[chargeIdx][e.getDC_sector()-1]->Fill( e.get3Momentum().Mag(), (e.getEpcal() + e.getEecin() + e.getEecout())/e.get3Momentum().Mag() );

				hBeta_p_e[0]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[0]->Fill(  e.getBeta() );
				hBeta_p_pi[0][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[0][chargeIdx]->Fill( pi_dummy.getBeta() );

				if( !anal.applyElectronFiducials( e ) ){ continue; }
				hBeta_p_e[1]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[1]->Fill( e.getBeta() );
				hBeta_p_pi[1][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[1][chargeIdx]->Fill( pi_dummy.getBeta() );
					
				if( !anal.applyElectronPCAL( e ) ){continue;}
				hBeta_p_e[2]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[2]->Fill( e.getBeta() );
				hBeta_p_pi[2][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[2][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyElectronEDep( e )){continue;}
				hBeta_p_e[3]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[3]->Fill( e.getBeta() );
				hBeta_p_pi[3][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[3][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyElectronSF( e )){continue;}
				hBeta_p_e[4]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[4]->Fill( e.getBeta() );
				hBeta_p_pi[4][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[4][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyElectronCorrelation( e )){continue;}
				hBeta_p_e[5]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[5]->Fill( e.getBeta() );
				hBeta_p_pi[5][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[5][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyElectronVertex( e )){continue;}
				hBeta_p_e[6]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[6]->Fill( e.getBeta() );
				hBeta_p_pi[6][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[6][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyElectronDetectorCuts( e )){continue;}
				hBeta_p_e[7]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[7]->Fill( e.getBeta() );
				hBeta_p_pi[7][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[7][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyPionDetectorFiducials( pi_dummy )){ continue ; }
				hBeta_p_e[8]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[8]->Fill( e.getBeta() );
				hBeta_p_pi[8][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[8][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyPionDetectorChi2( pi_dummy ) ){ continue; }
				hBeta_p_e[9]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[9]->Fill( e.getBeta() );
				hBeta_p_pi[9][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[9][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				if( !anal.applyPionDetectorVertex( pi_dummy, e ) ){ continue; }
				hBeta_p_e[10]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				hBeta_e[10]->Fill( e.getBeta() );
				hBeta_p_pi[10][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				hBeta_pi[10][chargeIdx]->Fill( pi_dummy.getBeta() );
					
				hFid_e[chargeIdx][0][1]->Fill( e.getDC_x1(), e.getDC_y1() ); 
				hFid_pi[chargeIdx][0][1]->Fill( pi_dummy.getDC_x1(), pi_dummy.getDC_y1() );
				
				hFid_e[chargeIdx][1][1]->Fill( e.getDC_x2(), e.getDC_y2() ); 
				hFid_pi[chargeIdx][1][1]->Fill( pi_dummy.getDC_x2(), pi_dummy.getDC_y2() );
				
				hFid_e[chargeIdx][2][1]->Fill( e.getDC_x3(), e.getDC_y3() ); 
				hFid_pi[chargeIdx][2][1]->Fill( pi_dummy.getDC_x3(), pi_dummy.getDC_y3() );
	
			}

			
		}
	}

	std::cout<< "Finished File loop! \n";
	
	std::cout<<"Writing tree to file\n";
   	outputFile->cd();
	hPositives->Write();
	hNegatives->Write();
	
	hPos_b_p->Write();
	hNeg_b_p->Write();
	for( int i = 0; i < 2; i++ ){
		hPCAL_WV[i]->Write();
		hEdep[i]->Write();
		hSF_corr[i]->Write();
		hVz_e[i]->Write();

		hVz_pi[i]->Write();
		hChi2[i]->Write();
	
		for( int j = 0; j < 3; j++ ){
			hFid_e[i][j][0]->Write();
			hFid_e[i][j][1]->Write();
			
			hFid_pi[i][j][0]->Write();
			hFid_pi[i][j][1]->Write();
		}
		
		for( int j = 0; j < 6; j++ ){
			hSF[i][j]->Write();
		}
	}

	for( int i = 0; i < 11; i++ ){
		hBeta_p_e[i]->Write();
		hBeta_p_pi[i][0]->Write();
		hBeta_p_pi[i][1]->Write();
	
		hBeta_e[i]->Write();
		hBeta_pi[i][0]->Write();
		hBeta_pi[i][1]->Write();
	}
	outputFile->Close();
	
	std::cout<<"Done!\n";

	auto finish = std::chrono::high_resolution_clock::now();
    	std::chrono::duration<double> elapsed = finish - start;

    	std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...................... END OF MAIN.................................//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug ){

	double beam_helicity = p_event->getHelicity();

	bool helFlip = true;
	if      (runnum>=11093 && runnum<=11283)    helFlip = false; // falls, 10.4 GeV period only
	else if (runnum>=11323 && runnum<=11571)    helFlip = false; // winter

	if (helFlip) {
		beam_helicity = -1 * beam_helicity;
	}
	return beam_helicity;
}



int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne){
	// find leading electron as the one with highest energy

	double  leading_e_E = 0;
	int     leading_e_index = 0;


	for (int eIdx=0; eIdx < Ne; eIdx++) {
		double M, P;   

		M = rp[eIdx]->getCalcMass();
		P = rp[eIdx]->par()->getP();
	
		double E_e = sqrt(M*M + P*P);

		if (E_e > leading_e_E) {
			leading_e_index = eIdx;
			leading_e_E     = E_e;
		}
	}

	return leading_e_index;
}




