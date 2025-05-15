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
#include "TChain.h"
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
#include "analyzer.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST
#define RUN_PATH _DATA

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

using namespace cutVals;
using namespace constants;

int main( int argc, char** argv){

	if( argc < 1 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Output Directory]\n";
		return -1;
	}

	TString out_name = argv[1];
	
	double energy_list[3] = {10.2, 10.4, 10.6};

	TChain * file_rec = new TChain("ePi");
	
	for( double energy : energy_list ){

		TFile * file = new TFile(Form("/volatile/clas12/users/jphelan/SIDIS/data/final_skims/%0.1f/final_skim.root", energy));	
		TTree * ePi = (TTree*)file->Get("ePi");
		TString runList;
		if(energy == 10.6){ 
			runList = (TString) RUN_PATH + "/runLists/good_runs_10-6.txt";
		}
		
		else if(energy == 10.4){ 
			runList = (TString) RUN_PATH+"/runLists/good_runs_10-4.txt";
		}
		
		else{
			runList = (TString) RUN_PATH + "/runLists/good_runs_10-2.txt";
		}

		cout<<"Run List : "<<runList<<endl;

		std::ofstream txtFile;
		txtFile.open( out_name + Form("/piCounts_%0.1f.txt",energy) );
		txtFile<<"Run\tN+\tN-\n";

		std::ifstream stream;
		stream.open(runList);
		int runnum;
		std::string runNum;
	
		
		while(std::getline(stream, runNum)){
			cout<<"Doing run "<<runnum<<std::endl;
			runnum = stoi(runNum);
			int pip = ePi->GetEntries(Form("pi.Charge>0 && isGoodPion && runnum==%i",runnum));
			int pim = ePi->GetEntries(Form("pi.Charge<0 && isGoodPion && runnum==%i",runnum));
			cout<<"Pi+ : "<<pip<<" & Pi- : "<<pim<<std::endl;
			txtFile<<runnum<<"\t"<<pip<<"\t"<<pim<<std::endl;
		}

		txtFile.close();
	}
}
