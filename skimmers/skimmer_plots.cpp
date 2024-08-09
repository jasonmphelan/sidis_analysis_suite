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
//#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
//#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
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

	if( argc < 7 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Beam energy]\n";
		cerr << "	[Run Type] [Single File?] [Inclusive (0, 1)] [Output File Name (no extension)]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
       	double Ebeam = atof(argv[2]); // [GeV]
	int RunType = atoi(argv[3]);
	int inclusive =atoi( argv[5]);
	int singleFile =atoi( argv[4]);
	TString outFileName = argv[6]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	
	if( singleFile != 0 && singleFile != 1 ){
		cout<<"Invalid entry for number of output files : "<<singleFile<<std::endl;
		cout<<"Outputting single file\n";
	}

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	anal.setAnalyzerLevel(0);
	anal.loadCutValues(-1, Ebeam);
	
	reader runReader;
	runReader.setNumFiles( nFiles);
	runReader.setRunType( RunType );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
       	runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;
	
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
	TTree * outTree;
	
	outputFile = new TFile(outFileName + ".root", "RECREATE");
	
	TH2F * hPositives = new TH2F("hPositives", ";p;#beta", 500, 0, 5, 500, 0.9, 1.01);
	TH2F * hNegatives = new TH2F("hNegativs", ";p;#beta", 500, 0, 5, 500, 0.9, 1.01);


	TH2F * hElectronFid = new TH2F( "hElectron_FID", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	TH2F * hElectronWV = new TH2F( "hElectron_WV", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	TH2F * hElectronMinE = new TH2F( "hElectron_MinE", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	TH2F * hElectronSF = new TH2F( "hElectron_SF", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	TH2F * hElectronSF_Corr = new TH2F( "hElectron_SF_Corr", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	TH2F * hElectronVt = new TH2F( "hElectron_Vt", ";p_{e};#beta", 500, 0, 5, 250, 0.9, 1.01 );

	TH2F * hPionFID[2];
	TH2F * hPionChi2[2];
	TH2F * hPionVt[2];
	TString chargeType[2] = {"Pip", "Pim"};	

	for( int i = 0; i < 2; i++ ){
		hPionFID[i] = new TH2F( Form("h%s_FID", chargeType[i]), ";p_{#pi};#beta", 500, 0, 5, 250, 0.9, 1.01 );
		hPionChi2[i] = new TH2F( Form("h%s_Chi2", chargeType[i]), ";p_{#pi};#beta", 500, 0, 5, 250, 0.9, 1.01 );
		hPionVt[i] = new TH2F( Form("h%s_Vt", chargeType[i]), ";p_{#pi};#beta", 500, 0, 5, 250, 0.9, 1.01 );
	}


	


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
           		if( RunType > 0 && event%1000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
           		if( RunType == 0 && event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
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
			std::vector<region_part_ptr> allParts = c12.getDetParticles();	
	
			if( RunType == 0 ){
				for( auto particle : allParts ){
						
					if( particle->par()->getCharge() > 0 ){
						hPositives->Fill(particle->par()->getP(),particle->getBeta() );	
					}
					if( particle->par()->getCharge() < 0 ){
						hNegatives->Fill( particle->par()->getP(),particle->getBeta());	
					}

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
			//Find good electrons
			int e_idx = GetLeadingElectron(electrons, Ne);	
			e.setElectron( Ebeam, electrons[e_idx]);
			double beta = e.getBeta();
			double mom = e.get3Momentum().Mag();

			////////////////Pion analysis/////////////////
			if(!anal.applyElectronFiducials( e )){ continue; };
			hElectronFid->Fill( mom, beta );
			
			if(!anal.applyElectronPCAL( e )){ continue; }
			hElectronWV->Fill( mom, beta );

			
			if(!anal.applyElectronEDep( e )){ continue; }
			hElectronMinE->Fill( mom, beta );

			
			if(!anal.applyElectronSF( e )){ continue; }
			hElectronSF->Fill( mom, beta );

			
			if(!anal.applyElectronCorrelation( e )){ continue; }
			hElectronSF_Corr->Fill( mom, beta );
 
			
			if(!anal.applyElectronVertex( e )){ continue; }
			hElectronVt->Fill( mom, beta );

			
			pion pi_dummy;
			//Find good pions			
			for(int i = 0; i < Npi; i++){
				if( inclusive == 1 ){ continue; }
				pi_dummy.Clear();
				pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
		
				beta = pi_dummy.getBeta();
				mom = pi_dummy.get3Momentum().Mag();
				int chargeIdx = (int)( pi_dummy.getCharge() < 0 );
				if(!anal.applyPionDetectorFiducials( pi_dummy )){continue;}
				hPionFID[chargeIdx]->Fill( mom, beta );
				
				if(!anal.applyPionDetectorChi2( pi_dummy )){continue;}
				hPionChi2[chargeIdx]->Fill( mom, beta );
				
				if(!anal.applyPionDetectorVertex( pi_dummy, e )){continue;}
				hPionVt[chargeIdx]->Fill( mom, beta );
			
					
			}

		}
	}

	std::cout<< "Finished File loop! \n";
	
	std::cout<<"Writing tree to file\n";
   	outputFile->cd();
	hPositives->Write();
	hNegatives->Write();

	hElectronFid->Write(); 
	hElectronWV->Write();  
	hElectronMinE->Write();  
	hElectronSF->Write();  
	hElectronSF_Corr ->Write(); 
	hElectronVt->Write();  


	for( int i = 0; i < 2; i++ ){
		hPionFID[i]->Write(); 
		hPionChi2[i]->Write(); 
		hPionVt[i]->Write(); 
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



