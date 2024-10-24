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
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"

#define HIST_PATH _HIST
#define CORR_PATH _DATA

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;

int main( int argc, char** argv){

	if( argc < 5 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [output] [x] [q2] [z] [p]\n";
		return -1;
	}
	
	int corrType = atoi(argv[1]);
	double x = atof(argv[2]);
	double q2 = atof(argv[3]);
	double z = atof(argv[4]);
	double p = atof(argv[5]);
	
	correctionTools corrector(1);
//	corrector.loadParameters();
	corrector.loadHistograms();
	//corrector.testHists();
	corrector.setKinematics(x, q2, z, p);

	cout<<"Correction Value : "<<corrector.getCorrectionFactor( corrType, 0 )<<endl;
}



