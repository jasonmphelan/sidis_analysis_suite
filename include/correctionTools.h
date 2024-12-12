#ifndef __CORRECTIONTOOLS__
#define __CORRECTIONTOOLS__

#include <sstream>
#include <fstream>
#include <iostream>
#include "TString.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "constants.h"
#include "cut_values.h"
#define CORR_PATH _DATA

using namespace constants;
using namespace cutVals;

class correctionTools{
public:
	correctionTools(int _mode);
	~correctionTools();
	
	void setCorrectionLevel( int _mode ){ mode = _mode; }
	void loadParameters();
	void loadHistograms();

	void setWeightName( TString newName ){ weight_name = newName; };
	void setPi2kName( TString newName ){ pi_To_kaon_name = newName; };
	void setK2piName( TString newName ){ kaon_To_pi_name = newName; };

	void setKinematics( double x, double q, double z, double p ); 

	double getCorrectionFactor( int type, int charge );
	double getCorrectionError( int type, int charge );
//	void setFilePaths(int corr, TString path);
	void printFilePaths();
	void testHists(){ accCorrection[1]->Print("all");}
	void loadFits();

	void setCorrectionFile(TString name){ weight_name = name; }
	void loadNewEnergy( double energy );
private:
	int 	mode;
	TString weight_name = "corrections_10.2_fit.root";
	TString kaon_To_pi_name = "corrections_k2pi_fit.root";
	TString pi_To_kaon_name = "corrections_pi2k_fit.root";

	double kin[4];
	
	TFile * weightFile;
	TFile * pi2kFile;
	TFile * k2piFile;

	
	TH3F *	accCorrection[2];
	TH3F *	binMigration[2];
	TH3F *	k_to_pi_Correction[2][4];
	TH3F *	pi_to_k_Correction[2][4];
	//TH3F *	rhoCorrection;
	
	TF1 * weightFit[2][bins_xB][bins_Q2];
	TF1 * pi2kFit[2][bins_xB][bins_Q2][4];
	TF1 * k2piFit[2][bins_xB][bins_Q2][4];
	
	double weight_Parameters[2][4][4][4];
	double k_To_pi_Parameters[2][4][5][4][4];
	double pi_To_k_Parameters[2][4][5][4][4];

	void loadWeightParameters();
	void loadkTopiParameters();
	void loadpiTokParameters();

};

#endif
