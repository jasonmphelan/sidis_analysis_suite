// Erez O. C., Oct-6, 2021
#ifndef __ANALYZER__
#define __ANALYZER__

#include "csv_reader.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include "TLorentzVector.h"
#include "TString.h"
#include "cut_values.h"
#include "constants.h"
#include "electron.h"
#include "pion.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"
#include "TH3F.h"

using namespace cutVals;
using namespace constants;

class analyzer{
public:
	analyzer(int _mode, int _torusBending_);
	~analyzer();
	void setAnalyzerLevel( int _mode ){ mode = _mode; }
	void loadElectron( electron *iE ){ e = iE; } 
	void loadPion( pion *iPi ){ pi = iPi; } 

	double	Chi2PID_pion_upperBound (Double_t p, Double_t C);
	double 	Chi2PID_pion_lowerBound (Double_t p, Double_t C);

	bool 	applyElectronDetectorCuts( electron e );
	bool 	applyElectronKinematicCuts( electron e );

	bool 	applyPionKinematicCuts( pion pi );
	bool 	applyPionDetectorCuts( pion pi, electron e );

	void	loadCutValues (int torusBending, double EBeam); //  -1 for In-bending, +1 for Out-bending
	//void	loadAcceptanceMatching (int torusBending=1); //  -1 for In-bending, +1 for Out-bending
	//void	loadCorrections (TFile corrFileName); //  -1 for In-bending, +1 for Out-bending
	//void	getCorrections ( electron e, pion pi);
	bool	applyAcceptanceMatching( pion pi, int dim );

	void	printCutValues ();
	void	SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;}
	double	ComputeLightConeFraction ( TLorentzVector p );
	double	calcQStar ( TVector3 eP3, TVector3 piP3, double Ebeam );
	void loadMatchingFunctions( TString fileName );
	void loadMatchingFunctions();
	int acceptance_match_3d_cont( double phi_part, double theta, double p, TF1 * fitFuncs[6][3], int chargeIdx);
	int acceptance_match_3d( double phi_part, double theta, double p, int charge);
	//double	fillDetectorHistograms();

private:
	int	fdebug;
	int	torusBending; // -1 for In-bending, +1 for Out-bending
	int 	mode;


	TH3F *	accCorrection;
	TH3F *	binCorrection;
	TH3F *	kaonCorrection;
	TH3F *	rhoCorrection;


	TF1 * match3d[6][3];

	e_pid epid;
	DCfid_SIDIS dcfid;
	electron * e;
	pion * pi;
		
	

};

#endif
