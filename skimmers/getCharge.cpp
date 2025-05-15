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

	if( argc < 4 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Beam energy]\n";
		cerr << "	[Run Type] [Output File Name (no extension)]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
    double Ebeam = atof(argv[2]); // [GeV]
	int RunType = atoi(argv[3]);
	TString outFile_name = argv[4];
	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	if( RunType == 4 ){
		anal.setAnalyzerLevel(0);
	}	
	else{ anal.setAnalyzerLevel(RunType); }
	anal.loadCutValues(-1, Ebeam);
	anal.loadSamplingFractionParams();
	anal.printCuts();

	reader runReader;
	runReader.setNumFiles( nFiles);
	runReader.setRunType( RunType );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
    runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;
	
	// Set output variables
	
	double accCharge = 0;

	std::ofstream txtFile;
	std::cout<<"outfile : "<<outFile_name+".txt"<<std::endl;
	txtFile.open(outFile_name + ".txt");
	txtFile<<"Run\tCharge\n";
	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		//create the event reader
		clas12reader c12(files.GetFileName(i).Data());

		c12.scalerReader();
		double charge = c12.getRunBeamCharge();
		double runnum;
		while((c12.next()==true)){
			//evnum  = c12.runconfig()->getEvent();
			runnum = c12.runconfig()->getRun();
			break;
		}
		txtFile<<c12.runconfig()->getRun()<<"\t"<<charge<<std::endl;
		accCharge += charge;
		
	}

	txtFile<<"0\t"<<accCharge<<std::endl;

	txtFile.close();

	std::cout<<"Done!\n";
	std::cout<<"Accumulated Charge : "<<accCharge<<std::endl;
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



