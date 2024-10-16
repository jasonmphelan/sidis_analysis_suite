#ifndef __CORRECTIONTOOLS__
#define __CORRECTIONTOOLS__

#include <sstream>
#include <fstream>
#include <iostream>
#include "TString.h"
#include "TH3F.h"
#include "TFile.h"
#define CORR_PATH _DATA


class correctionTools{
public:
	correctionTools(int _mode);
	~correctionTools();
	
	void setCorrectionLevel( int _mode ){ mode = _mode; }
	void loadParameters();
	void loadHistograms();

	void setKinematics( double x, double q, double z, double p ); 

	double getCorrectionFactor( int type, int charge );
	double getCorrectionError( int type, int charge );
	void setFilePaths(int corr, TString path);
	void printFilePaths();
	void testHists(){ accCorrection[1]->Print("all");}
private:
	int 	mode;
	TString weight_name = "corrections.root";
	TString kaon_To_pi_name = "corrections_k_To_pi.root";
	TString pi_To_kaon_name = "corrections_pi_To_k.root";

	double kin[4];
	
	TFile * weightFile;
	TFile * pi2kFile;
	TFile * k2piFile;

	
	TH3F *	accCorrection[2];
	TH3F *	binMigration[2];
	TH3F *	k_to_pi_Correction[2][4];
	TH3F *	pi_to_k_Correction[2][4];
	//TH3F *	rhoCorrection;
	

	double weight_Parameters[2][4][4][4];
	double k_To_pi_Parameters[2][4][5][4][4];
	double pi_To_k_Parameters[2][4][5][4][4];

	void loadWeightParameters();
	void loadkTopiParameters();
	void loadpiTokParameters();

};

#endif
