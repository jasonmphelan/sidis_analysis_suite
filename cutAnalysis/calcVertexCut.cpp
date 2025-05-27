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
#include "TArray.h"

using namespace clas12;
using namespace constants;

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
Cuts applied in order
1) PCAL FID
2) PCAL MIN
3) SF
4) SF Correlation
5) Vz_e
6) Electron DC Fiducials
7) Vz_pi
8) Chi2
9) DC Fiducials
*/

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [outFile name (no extension)] [Beam energy (0 for all energies)]\n";
		return -1;
	}
	
       	double Ebeam = atof(argv[2]); // [GeV]
	TString outFile_name = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	TFile * outFile = new TFile(outFile_name, "RECREATE");
	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 && Ebeam != 0 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	std::cout<<"Loaded cut values\n";
	anal.setAnalyzerLevel(1);
	anal.loadCutValues(-1, Ebeam);
	
	reader runReader;
	runReader.setNumFiles( 1 );
	runReader.setRunType( 0 );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
       	runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;

	//Declare hists
	TString spec[2] = {"pip", "pim"};

	TH1F * hVz_e[2];

	TH1F * hVz_pi[2];
	
	for( int i = 0; i < 2; i++ ){
		hVz_e[i] = new TH1F( "hVz_e_"+spec[i], ";V_{Z}^{e} [cm];Counts [a.u.]", 400, -30, 30);

		hVz_pi[i] = new TH1F( "hVz_"+spec[i], ";V_{Z}^{#pi} - V_{Z}^{e} [cm];Counts [a.u.]", 400, -30, 30);
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
	//TFile * outputFile = new TFile(outFile_name + ".root", "RECREATE");
	int nFiles = 0;
	int RunType = 0;
	int inclusive = 0;

	//Set event count array
	double counts[2][20] = {0}; //[charge][cut level]

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
			
			// Get Particles By PID
			electrons   = c12.getByID( -11   );
			pipluses    = c12.getByID( 211  );
			piminuses   = c12.getByID(-211  );
			
			pions = pipluses;
			pions.insert( pions.end(), piminuses.begin(), piminuses.end() );

			Ne      = electrons.size();
			Npi 	= pions.size();
			
			if( Ne < 1 ){ continue; } //Keep only events with one electron...
			
			//////////////electron analysis////////////////////
			int e_idx = GetLeadingElectron(electrons, Ne);	
			e.setElectron( Ebeam, electrons[e_idx]);
			////////////////Pion analysis/////////////////
			
			//Find good pions			
			//if( !anal.applyElectronFiducials( e ) ){ continue; }
			//if( !anal.applyElectronPCAL( e ) ){continue;}
			//if( !anal.applyElectronEDep( e )){continue;}
			//if( !anal.applyElectronSF( e )){continue;}
			//if( !anal.applyElectronCorrelation( e )){continue;}
			hVz_e[0]->Fill( e.getVt().Z() );
			
			//if( !anal.applyElectronVertex( e )){continue;}	
			if( Npi == 0 && inclusive != 1 ){ continue; }	
			pion pi_dummy;

			for(int i = 0; i < Npi; i++){

				pi_dummy.Clear();
				pi_dummy.setMCPion( (bool)RunType );
				pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
	
				int chargeIdx = (int)(pi_dummy.getCharge() < 0);

				//if( !anal.applyPionDetectorFiducials( pi_dummy )){ continue ; }

				hVz_pi[chargeIdx]->Fill( pi_dummy.getVt().Z() - e.getVt().Z() );
				//if( !anal.applyPionDetectorVertex( pi_dummy, e ) ){ continue; }
				
				//counts[chargeIdx][10]++;
			}

			
		}
	}

	std::cout<< "Finished File loop! \n";
	
	std::vector<TH1F *> hist_list;
	hist_list.push_back(hVz_e[0]); 
	hist_list.push_back(hVz_pi[0]); 
	hist_list.push_back(hVz_pi[1]); 
	outFile->cd();
	for( TH1F * h : hist_list ){
	
	

		//Determine cut values
		double cutMin = -999;
		double cutMax = 999;
		bool aboveMax = false;
		bool belowMax = true;


		double maxi = h->GetMaximum();
		int i = 1;
		
		while( !aboveMax ){
			double ev = h->GetBinContent(i);
			if( ev > 0.25*maxi ){
				cutMin = h->GetBinCenter(i);
				belowMax = false;
				aboveMax = true;
			}
			i++;
		}
		while( !belowMax ){
			double ev = h->GetBinContent(i);
			if( ev < 0.25*maxi ){
				cutMax = h->GetBinCenter(i);
				belowMax = true;
				aboveMax = false;
			}
			i++;
		}
			
		cout<<"Min e vertex : "<<cutMin<<std::endl;
		cout<<"Max e vertex : "<<cutMax<<std::endl;
		
		
		h->Write();
	}

	outFile->Close();
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




