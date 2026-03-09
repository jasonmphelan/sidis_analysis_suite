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
	void setWeightName4D( TString newName ){ weight_name_4D = newName; };
	void setPi2kName( TString newName ){ pi_To_kaon_name = newName; };
	void setK2piName( TString newName ){ kaon_To_pi_name = newName; };

	void setKinematics( double x, double q, double z, double p ); 
	void set4dBin(int inBin){ bin_4d = inBin; }
	void setN4dBins(int nBins ) { nBins_4d = nBins; }
	double getCorrectionFactor( int type, int charge );
	double getCorrectionError( int type, int charge );
//	void setFilePaths(int corr, TString path);
	void printFilePaths();
	void testHists(){ pi_to_k_Correction[0][2]->Print("all");}
	void loadFits();

	void setCorrectionFile(TString name){ weight_name = name; }
	void loadNewEnergy( double energy );
	// Probe 4D correction file: returns bin count and populates corr_var_min/max/name
	int probeN4dBins();
	double getCorrVarMin()  { return corr_var_min; }
	double getCorrVarMax()  { return corr_var_max; }
	TString getCorrVarName(){ return corr_var_name; }

	double getX(){ return kin[0];}
	double getQ2(){ return kin[1];}
	double getZ(){ return kin[2];}
	double getP(){ return kin[3];}



private:
	int 	mode;
	TString weight_name = "corrections_10.2_AN_test.root";
	TString weight_name_4D = "corrections_10.2_phi_q.root";
	TString kaon_To_pi_name = "corrections_k2pi_AN.root";
	TString pi_To_kaon_name = "corrections_pi2k_AN.root";

	//TFile * weightFile = new TFile((TString) _DATA + "/correctionFiles/"+ weight_name);
	//TFile * pi2kFile = new TFile((TString) _DATA + "/correctionFiles/"+ pi_To_kaon_name);
	//TFile * k2piFile= new TFile((TString) _DATA + "/correctionFiles/"+ kaon_To_pi_name);

	double kin[4];
	int bin_4d = -1;
	int nBins_4d = 0;
	double corr_var_min  = 0;
	double corr_var_max  = 1;
	TString corr_var_name = "";

	TH3F *	accCorrection[2];
	TH3F *	binMigration[2];
	TH3F *	mcCorrection[2];
	
	TH3F *	mcCorrection4D[2][10];

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
	void loadWeightParameters_4D();
	void loadkTopiParameters();
	void loadpiTokParameters();

};

#endif
