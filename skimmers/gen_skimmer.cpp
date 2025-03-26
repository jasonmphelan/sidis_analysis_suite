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
//#include "electron.h"
//#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
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

	if( argc < 6 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Beam energy]\n";
		cerr << "	[Single File?] [Inclusive (0, 1)] [Output File Name (no extension)]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
       	double Ebeam = atof(argv[2]); // [GeV]
	int singleFile =atoi( argv[3]);
	int inclusive =atoi( argv[4]);
	TString outFileName = argv[5]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
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
	anal.setAnalyzerLevel(2);
	anal.loadCutValues(-1, Ebeam);
	
	reader runReader;
	runReader.setNumFiles( nFiles);
	runReader.setRunType( 1 );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
       	runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;
	
	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );

	// Electron Variables
	genElectron e_gen;

	// Pion Variables
	std::vector<genPion> pi_gen ;
	
	std::vector<int> electronsMC;
	std::vector<int> pionsMC, pipsMC, pimsMC;
	
	// Set Output file and tree
	TFile * outputFile;
	TTree * outTree;
	
	if( singleFile == 1 ){
		outputFile = new TFile(outFileName + ".root", "RECREATE");
		outTree = new TTree("ePi", "(e,e'pi) event  information");

		// Set output branches
		cout<<"Declare trees"<<endl;
		outTree->Branch("runnum", &runnum);
		outTree->Branch("torus", &torusBending);
		outTree->Branch("evnum", &evnum);
		outTree->Branch("Ebeam", &Ebeam);
		outTree->Branch("beam", &beam);
		
		outTree->Branch("Ne", &Ne);
			
		outTree->Branch("Npi", &Npi);
		outTree->Branch("Npips", &Npips);
		outTree->Branch("Npims", &Npims);

		outTree->Branch("e_gen", &e_gen);
		outTree->Branch("pi_gen", &pi_gen);
	}
	
	double accCharge = 0;
	int goodElectron = 0;

	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		if( singleFile == 0 ){
			cout<<"Declare output file\n";	
			outputFile = new TFile(outFileName + Form("_%i.root", runReader.getRunNum(i) ), "RECREATE");
			outTree = new TTree("ePi", "(e,e'pi) event  information");

			// Set output branches
			cout<<"Declare tree"<<endl;
			outTree->Branch("runnum", &runnum);
			outTree->Branch("torus", &torusBending);
			outTree->Branch("evnum", &evnum);
			outTree->Branch("Ebeam", &Ebeam);
			outTree->Branch("beam", &beam);
		
			outTree->Branch("Ne", &Ne);
		
			outTree->Branch("Npi", &Npi);
			outTree->Branch("Npips", &Npips);
			outTree->Branch("Npims", &Npims);

			outTree->Branch("e_gen", &e_gen);
			outTree->Branch("pi_gen", &pi_gen);
		}

		//create the event reader
		clas12reader c12(files.GetFileName(i).Data());
		auto mcparts = c12.mcparts();		
		int NeventsTotal = c12.getReader().getEntries();       
    	   	int event = 0;	

    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( event%1000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			//if(event > 10000){break;}	
			//Get run and event info	
			
			evnum  = c12.runconfig()->getEvent();
			runnum = c12.runconfig()->getRun();
			///////////////////////////Initialize variables//////////////////////////////////////////////	
			//electronsMC.clear();
			pionsMC.clear();
			pipsMC.clear();
			pimsMC.clear();
			
			e_gen.Clear();
			pi_gen.clear();

			Ne = Npi = Npips = Npims = 0;

			/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
			
			// Get Particles By PID
			
			
			int nMcPart = c12.mcevent()->getNpart();	
			int mcId;
			for ( int i = 1; i<nMcPart; i++ ){
				mcId = c12.mcparts()->getPid(i);
				switch (mcId){
					case 11:
						electronsMC.push_back(i);
						break;
					case 211:
						pipsMC.push_back(i);
						break;
					case -211:
						pimsMC.push_back(i);
						break;
				}
			}
			pionsMC = pipsMC;
			pionsMC.insert( pionsMC.end(), pimsMC.begin(), pimsMC.end() );
			Npi 	= pionsMC.size();
			Ne      = electronsMC.size();
			
			//if( Ne < 1 ){ continue; } //Keep only events with one electron...
			if( Npi == 0 && inclusive != 1 ){ continue; }	
			
           		//if( event%1000 == 0){cout<<"Found Particles for Event: "<<event<< "/"<<NeventsTotal<<endl; }
			
			//////////////electron analysis////////////////////
			//Find good electrons
			int e_idx = electronsMC[0];//no ambiguity with MC electron
			mcparts->setEntry(e_idx);
			e_gen.setMomentum( mcparts );
			e_gen.setKinematicInformation( Ebeam, mcparts );
			if( !anal.applyElectronKinematicCuts( e_gen )){continue;}
		
           		//if( event%1000 == 0){cout<<"Analyzed e- for Event: "<<event<< "/"<<NeventsTotal<<endl; }
			////////////////Pion analysis/////////////////
			
			genPion genPi_dummy;
			//Find good pions			
			for(int i = 0; i < Npi; i++){
				if( inclusive == 1 ){ continue; }
				genPi_dummy.Clear();
				mcparts->setEntry( pionsMC[i] );
				genPi_dummy.setMomentum( mcparts );
				genPi_dummy.setKinematicInformation( e_gen.getQ(),e_gen.get4Momentum(), mcparts );
				if( !anal.applyPionKinematicCuts( genPi_dummy ) ) {continue;}
				
				pi_gen.push_back(genPi_dummy);
					
			}

			//goodElectron++;
			if(pi_gen.size() == 0 && inclusive == 0){continue;}	
			outTree->Fill();
			
		}
		std::cout<<"Finished File!\n";
   		if( singleFile == 0 ){
			std::cout<<"Writing file!\n";
			outputFile->cd();
    			outTree->Write();
			outputFile->Close();
		}
	}

	std::cout<< "Finished File loop! \n";
	
	if( singleFile == 1 ){
		std::cout<<"Writing tree to file\n";
   		outputFile->cd();
    		outTree->Write();
		outputFile->Close();
	}
	
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



