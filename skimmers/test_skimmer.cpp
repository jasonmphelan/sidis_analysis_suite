/*
 * THINGS TO DO:
 * 1) Write gen-reco matching
 * 2) Figure out Q2, W issue 
*/

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
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

using namespace clas12;
using namespace constants;

void getRunFiles(TString path, TString runList, clas12root::HipoChain &files, int runType, int nFiles, int beamType);

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader * c12);

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList );

int FindMatch(TLorentzVector p, mcpar_ptr mcparts, std::vector<int> part_list, int PID);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 7 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Debug mode (0, 1)] [Beam energy]\n";
		cerr << "	[Run Type] [Inclusive (0, 1)] [Output File Name]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
        int fdebug = atoi(argv[2]);
       	double Ebeam = atof(argv[3]); // [GeV]
	int RunType = atoi(argv[4]);
	int inclusive =atoi( argv[5]);
	TString outFileName = argv[6]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	
	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	anal.setAnalyzerLevel(0);
	anal.loadCutValues(-1, Ebeam);

	
	// Set input file list    
	// Only seems to work with hipo chain method...

	TString DataPath;
	string runList;
	clas12root::HipoChain files;


	setDataPaths( RunType, Ebeam, DataPath, runList );	
	getRunFiles( DataPath, runList, files, RunType, nFiles, (int)( ( Ebeam - 10.2 )/0.2 ) );

	cout<<"Set output files"<<endl;
	// Set Output file and tree
	TFile * outputFile = new TFile(outFileName, "RECREATE");
	TTree * outTree = new TTree("ePi", "(e,e'pi) event  information");


	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );

	// Electron Variables
	electron e;
	
	// Pion Variables
	std::vector<pion> pi ;
	
	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted
	
	// Set output branches
	cout<<"Declare trees"<<endl;
	outTree->Branch("runnum", &runnum);
	outTree->Branch("torus", &torusBending);
	outTree->Branch("evnum", &evnum);
	outTree->Branch("Ebeam", &Ebeam);
	outTree->Branch("beam", &beam);
	
	outTree->Branch("Ne", &Ne);
	outTree->Branch("e", &e);
		
	outTree->Branch("Npi", &Npi);
	outTree->Branch("Npips", &Npips);
	outTree->Branch("Npims", &Npims);

	outTree->Branch("pi", &pi);
	
	double accCharge = 0;
	int goodElectron = 0;


	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	if (fdebug) std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;
	   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		//create the event reader
		//clas12reader c12(files->At(i)->GetTitle(),{0});
		clas12reader c12(files.GetFileName(i).Data());
		auto mcparts = c12.mcparts();		

		int NeventsTotal = c12.getReader().getEntries();       
    	   	int event = 0;
			
    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( RunType > 0 && event%1000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
           		if( RunType == 0 && event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			//if(event > 10000){break;}	
			//Get run and event info	
			
			runnum = c12.runconfig()->getRun();
			evnum  = c12.runconfig()->getEvent();
			
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
			if( !anal.applyElectronDetectorCuts( e )){continue;}
		
			////////////////Pion analysis/////////////////
			
			pion pi_dummy;
			//Find good pions			
			for(int i = 0; i < Npi; i++){
				if( inclusive == 1 ){ continue; }
				pi_dummy.Clear();
			
				pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
				if( !anal.applyPionDetectorCuts( pi_dummy, e ) ) {continue;}
				pi.push_back(pi_dummy);

			}

			//goodElectron++;
			if(pi.size() == 0 && inclusive == 0){continue;}
			
	
			outTree->Fill();
			
		}
	
	}
       
   	outputFile->cd();
    	outTree->Write();
	outputFile->Close();


	auto finish = std::chrono::high_resolution_clock::now();
    	std::chrono::duration<double> elapsed = finish - start;

    	std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...................... END OF MAIN.................................//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList ){

	if(runType == 1){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/job_%i/output/%i-",eBeam, job, job);//reco_out/clasdis_rec_", eBeam);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    		//dataPath = "/volatile/clas12/osg/jphelan/job_7393/output/7393-";//reco_out/clasdis_rec_", eBeam);
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/reco_out/clasdis_rec_", eBeam);
    		//dataPath = "/work/clas12/users/jphelan/SIDIS_at_BAND/generator/dst_";
	}
	
	else if( runType == 2 ){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/claspyth/%.1f/proton/reco_out/claspyth_rec_", eBeam);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    	}
	
	else if( runType == 3 ){
		dataPath = "/volatile/clas12/users/jphelan/SIDIS/data/rho_skims/reco_out/reco_rho_";
    	}
	

	else{ //runType == 0 is data
		TString path_temp;
	
		if(eBeam == 10.6){ 
			runList = "/work/clas12/users/jphelan/SIDIS_at_BAND/macros/runlists/good_runs_10-6.txt";
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else if(eBeam == 10.4){ 
			runList = "/work/clas12/users/jphelan/SIDIS_at_BAND/macros/runlists/good_runs_10-4.txt";
			path_temp = "spring2020/torus-1/pass2/v1/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else{
			runList = "/work/clas12/users/jphelan/SIDIS_at_BAND/macros/runlists/good_runs_10-2-final.txt";
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
    		}
    		
		dataPath = "/cache/clas12/rg-b/production/recon/"+path_temp;	
    	}
}

void getRunFiles(TString path, TString runList, clas12root::HipoChain &files, int runType, int nFiles, int beamType){
	std::ifstream stream;
	stream.open(runList);
	string runNum;
	TString inFile;
	if(runType == 0){
		int i = 0;
		while(std::getline(stream, runNum)){
			if(nFiles != 0 && i >= nFiles ) break;
			TString run(runNum);
			inFile = path + run+".hipo";
			cout<<inFile<<std::endl;
			files.Add(inFile.Data());		
			i++;
		}
	}
	else if (runType == 1){
		for( int j = 0; j < nRuns[beamType]; j++ ){
			for( int i = 0; i < 75000; i++){
				if( nFiles != 0 && i >= nFiles ) break;
				inFile = path + Form("%i/output/%i-%i.hipo", monteCarloRuns[beamType][j], monteCarloRuns[beamType][j], i+1);
				if( gSystem->AccessPathName(inFile) ) continue;
				files.Add(inFile.Data());
			}
		}
	}
	else {
		for( int i = 0; i < 18; i++){
			if( nFiles != 0 && i >= nFiles ) break;
			inFile = path + Form("%i.hipo", i);
			//if( gSystem->AccessPathName(inFile) ) continue;
			files.Add(inFile.Data());
		}
	}

	//else{
	//	for( int i = 0; i < 1000; i++){
	//		if( nFiles != 0 && i >= nFiles ) break;
	//		inFile = path + Form("%i.hipo", i+1);
	//		if( gSystem->AccessPathName(inFile) ) continue;
	//		files.Add(inFile.Data());
	//	}
	//}
}





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

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader * c12){
	// find leading electron as the one with highest energy
	// generator level
	// note overloaded function
	
	double  leading_e_E = 0;
	int     leading_e_index = 0;


	for (int eIdx=0; eIdx < Ne; eIdx++) {
		double M, P2;   

		//For generator level, need mc particle info
		M = c12->mcparts()->getMass(eIdx);
		P2 = c12->mcparts()->getPx(eIdx)*c12->mcparts()->getPx(eIdx) + 
		c12->mcparts()->getPy(eIdx)*c12->mcparts()->getPy(eIdx) +
		c12->mcparts()->getPz(eIdx)*c12->mcparts()->getPz(eIdx);
	
		double E_e = sqrt(M*M + P2);

		if (E_e > leading_e_E) {
			leading_e_index = eIdx;
			leading_e_E     = E_e;
		}
	}

	return leading_e_index;
}



int FindMatch(TLorentzVector p, mcpar_ptr mcparts, std::vector<int> part_list, int PID){

	double temp_min_dp = 9999;
	double p_idx = -1;
	double temp_min_dtheta = 9999;
	double theta_idx = -1;

	for ( int i : part_list ){
		mcparts->setEntry(i);

		if( mcparts->getPid() != PID ){ continue; }		

		if( abs( p.Vect().Mag() - mcparts->getP() ) < temp_min_dp ){	
			temp_min_dp = abs( p.Vect().Mag() - mcparts->getP() );
			p_idx = i;
		}		
		if( abs( p.Theta() - mcparts->getTheta() ) < temp_min_dtheta ){	
			temp_min_dtheta = abs( p.Theta() - mcparts->getTheta() );
			theta_idx = i;
		}	
	}
	
	if( theta_idx != p_idx ){
		return -1;
	}
	
	return p_idx;
}

