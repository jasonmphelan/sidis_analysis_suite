#include "reader.h"
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "TF1.h"
#include "TString.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "HipoChain.h"
#include <fstream>
#include <unistd.h>
#include "constants.h"


#define RUN_PATH _DATA

using namespace std;
using namespace constants;

reader::reader() {}

reader::~reader() {}


int reader::getRunNum( int i ){
	if( runType != 0 ){return i;}
	else{ return runNums[i]; }
	return 0;
}
		


void reader::readRunFiles(clas12root::HipoChain &fileList){
	if( EBeam == 0 ){
		readRunFilesAllE( fileList );
	}

	else{
		setDataPaths();
		getRunList();
		getRunFiles(fileList);
	}
}

void reader::readRunFilesAllE(clas12root::HipoChain &fileList){
	setEnergy(10.2);
	readRunFiles(fileList);
	setEnergy(10.4);
	readRunFiles(fileList);
	setEnergy(10.6);
	readRunFiles(fileList);

}
	

void reader::setDataPaths(){

	if(runType == 1 || runType == 2){
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
	}
	

	else{ //runType == 0 is data
		TString path_temp;
	
		if(EBeam == 10.6){ 
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else if(EBeam == 10.4){ 
			path_temp = "spring2020/torus-1/pass2/v1/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else{
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
    		}
    		
		dataPath = "/cache/clas12/rg-b/production/recon/"+path_temp;	
    	}
}

void reader::getRunList(){
	if(EBeam == 10.6){ 
		runList = (TString) RUN_PATH + "/runLists/good_runs_10-6.txt";
	}
	
	else if(EBeam == 10.4){ 
		runList = (TString) RUN_PATH+"/runLists/good_runs_10-4.txt";
	}
	
	else{
		runList = (TString) RUN_PATH + "/runLists/good_runs_10-2.txt";
    	}

	cout<<"Run List : "<<runList<<endl;
}


void reader::getRunFiles( clas12root::HipoChain &files){

	//clas12root::HipoChain files;
	std::ifstream stream;
	stream.open(runList);
	string runNum;
	TString inFile;
	int beamType =	(int) ( (EBeam - 10.2)/.2 );
	if(!stream ){ cout<<"fAILE TO open runlist\n";}
	if(runType == 0){
		int i = 0;
		while(std::getline(stream, runNum)){
			if(nFiles != 0 && i >= nFiles ) break;
			TString run(runNum);
			inFile = dataPath + run+".hipo";
			cout<<inFile<<std::endl;
			runNums.push_back( atoi( run ) );
			files.Add(inFile.Data());		
			i++;
		}
	}
	else if (runType == 1){
		for( int j = 0; j < nRuns[beamType]; j++ ){
			for( int i = 0; i < 75000; i++){
				if( nFiles != 0 && i >= nFiles ) break;
				inFile = dataPath + Form("%i/output/%i-%i.hipo", monteCarloRuns[beamType][j], monteCarloRuns[beamType][j], i+1);
				if( gSystem->AccessPathName(inFile) ) continue;
				files.Add(inFile.Data());
			}
		}
	}
	else {
		for( int i = 0; i < 18; i++){
			if( nFiles != 0 && i >= nFiles ) break;
			inFile = dataPath + Form("%i.hipo", i);
			if( gSystem->AccessPathName(inFile) ) continue;
			files.Add(inFile.Data());
		}
	}


	//return files;
}


void reader::getSkimsByName( TChain * chain, TString name ){

	std::ifstream stream;
	stream.open(runList);
	string runNum;
	TString inFile;
	int beamType =	(int) ( (EBeam - 10.2)/.2 );
	if(!stream ){ cout<<"fAILE TO open runlist\n";}
	if(runType == 0){
		int i = 0;
		while(std::getline(stream, runNum)){
			if(nFiles != 0 && i >= nFiles ) break;
			TString run(runNum);
			inFile = name + Form("_%i", atoi(run)) + ".root";
			cout<<"Adding : "<<inFile<<endl;
			if( gSystem->AccessPathName(inFile) ) continue;
			chain->Add(inFile);		
			i++;
		}
	}
	else if (runType == 1){
		//for( int j = 0; j < nRuns[beamType]; j++ ){
			for( int i = 0;  i < 75000; i++){
				if( nFiles != 0 && i >= nFiles ) break;
				inFile = name + ".root";
				cout<<"File "<<i<<" : "<<inFile<<endl;
				if( gSystem->AccessPathName(inFile) ) continue;
				chain->Add(inFile);
			}
		//}
	}
	else{ cout<<"(Currently) invalid run type... no files added\n"; }
//	else {
//		for( int i = 0; i < 18; i++){
//			if( nFiles != 0 && i >= nFiles ) break;
//			inFile = dataPath + Form("%i.hipo", i);
//			if( gSystem->AccessPathName(inFile) ) continue;
//			files.Add(inFile.Data());
//		}
//	}


	//return files;
}
void reader::getRunSkimsByName( TChain * chain, TString name ){
	if( EBeam == 0 ){
		getRunSkimsAllEnergy(chain, name);
	}
	else{
		getRunList();
		getSkimsByName( chain, name );
	}
}
void reader::getRunSkimsAllEnergy( TChain * chain, TString name ){
	setEnergy(10.2);
	getRunSkimsByName(chain, name + Form("%i", 10.2) + "run_skim");
	setEnergy(10.4);
	getRunSkimsByName(chain, name+ Form("%i", 10.2) + "run_skim");
	setEnergy(10.6);
	getRunSkimsByName(chain, name+ Form("%i", 10.2) + "run_skim");
}
