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
#include "clas12databases.h"
//#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
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

	if( argc < 2 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output File Name (no extension)]\n";
		return -1;
	}
	
	int nFiles = 0;//atoi(argv[1]); //set 0 to loop over all files,
	TString outFileName = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	TH1F * hCharge_run[3];

	double runMin[3] = {6400, 6100, 11200};
	double runMax[3] = {6600, 6500, 11600};

	for ( int i = 0; i < 3; i++ ){
		hCharge_run[i] = new TH1F(Form("hCharge_run_%.1f", 10.2 + i*.2), "", (int)(runMax[i] - runMin[i]), runMin[i], runMax[i] );
	}
	std::vector<double> energy;
	energy.push_back(10.2);
	energy.push_back(10.4);
	energy.push_back(10.6);

	TFile * outputFile = new TFile(outFileName + ".root", "RECREATE");
	outputFile->cd();
	for( double Ebeam : energy ){

		// Read cut values
		double torusBending = -1; //outBending = -1, inBending = 1
		analyzer anal(0, torusBending);
		anal.setAnalyzerLevel(0);
		anal.loadCutValues(-1, Ebeam);
		
		reader runReader;
		runReader.setNumFiles( nFiles);
		runReader.setRunType( 0 );
		runReader.setEnergy( Ebeam );
		
		clas12root::HipoChain files;
		runReader.readRunFiles(files);

		cout<<"Set output files"<<endl;
		
		int eIdx = ((int) ( (Ebeam - 10.2)/.2 ) );
		
		

		////////////////////////////////////Begin file loop////////////////////////////////////////////////////
		for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
			std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;   	
			//Only skim desired number of files
			if(nFiles != 0 && i > nFiles){break;}	
			int evNum = -1;
			int runNum = 0;
			double runCharge = 0;


			clas12reader c12(files.GetFileName(i).Data());

			c12.scalerReader();
			runCharge = c12.getRunBeamCharge();



			while((c12.next()==true)){
				evNum++;
				
				//create the event reader

				if(evNum){
					runNum = (int) c12.runconfig()->getRun();
			
					cout<<"For run number : "<<runNum;
					cout<<", acc charge = "<< runCharge<<std::endl;


					hCharge_run[eIdx]->SetBinContent( runNum - runMin[eIdx] + 1, runCharge ); 
					break;
				}
			}
		}

		hCharge_run[eIdx]->Write();

		std::cout<< "Finished File loop! \n";
		
		
		std::cout<<"Done!\n";

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;

		std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
	}
	outputFile->Close();
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



