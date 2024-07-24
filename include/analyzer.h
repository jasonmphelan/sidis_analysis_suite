// Erez O. C., Oct-6, 2021
#ifndef __SIDISatBAND_auxiliary_H__
#define __SIDISatBAND_auxiliary_H__

#include "csv_reader.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include "TLorentzVector.h"
#include "TString.h"
#include "cutValues.h"
#include "electrons.h"
#include "pions.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"

using namespace cutValues;
using namespace constants;

class analyzer {
public:
	analyzer(int _mode, int _torusBending_=-1);
	~analyzer();

	void loadElectron( electron *iE ){ e = iE; } 
	void loadPion( pion *iPi ){ pi = iPi; } 

	double	Chi2PID_pion_upperBound (Double_t p, Double_t C);
	double 	Chi2PID_pion_lowerBound (Double_t p, Double_t C);

	bool 	applyElectronDetectorCuts( electron e );
	//bool 	applyElectronKinematics( electron e );

	bool 	applyPionKinematics( pion pi );
	//bool 	applyPionDetectorCuts( pion pi );

	void	loadCutValues (int torusBending=1); //  -1 for In-bending, +1 for Out-bending
	//void	loadAcceptanceMatching (int torusBending=1); //  -1 for In-bending, +1 for Out-bending
	//void	loadCorrections (TFile corrFileName); //  -1 for In-bending, +1 for Out-bending
	//void	getCorrections ( electron e, pion pi);
	//void	applyAcceptanceMatching( pion pi );	
	void	printCutValues ();
	void	SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;};
	//double	ComputeLightConeFraction ( TLorentzVector p );
	//double	calcQStar ( TVector3 eP3, TVector3 piP3, double Ebeam );
	//double	fillDetectorHistograms();

private:
	int	fdebug;
	int	torusBending; // -1 for In-bending, +1 for Out-bending
	
	TH3F *	accCorrection;
	TH3F *	binCorrection;
	TH3F *	kaonCorrection;
	TH3F *	rhoCorrection;
	
	e_pid ePID;
	DCfid_SIDIS DC_fid;
	electron * e;
	pion * pi;
		
	

};

#endif
