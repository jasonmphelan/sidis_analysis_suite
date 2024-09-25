#ifndef __CORRECTIONTOOLS__
#define __CORRECTIONTOOLS__

#include <sstream>
#include <fstream>
#include <iostream>
#include "TString.h"
#include "TH3F.h"

#define CORR_PATH _DATA

using namespace cutVals;
using namespace constants;

class correctionTools{
public:
	correctionTools(int _mode);
	~correctionTools();
	
	void setCorrectionLevel( int _mode ){ mode = _mode; }
	void loadParameters();

	void setKinematics( double x, double q, double z, double p ); 

	double getCorrectionFactor( int type, int charge );
	void setFilePaths(int corr, TString path);

private:
	int 	mode;
	TString weight_name = "corrections.root";
	TString kaon_To_pi_name = "corrections_kaon_To_pi.root";
	TString pi_To_kaon_name = "corrections_pi_To_kaon.root";

	double kin[4];

	TH3F *	accCorrection;
	TH3F *	binCorrection;
	TH3F *	k_to_pi_Correction;
	TH3F *	pi_to_k_Correction;
	TH3F *	rhoCorrection;

	double weight_Params[2][4][4][4];
	double k_to_pi_Params[2][4][4][4][4];
	double pi_to_k_Params[2][4][4][4][4];

};

#endif
