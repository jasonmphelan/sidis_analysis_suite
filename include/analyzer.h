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
#include "genElectron.h"
#include "pion.h"
#include "genPion.h"
#include "e_pid.h"
#include "DCfid_SIDIS.h"
#include "TH3F.h"
#include "mcparticle.h"

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
	bool 	applyElectronKinematicCuts( genElectron e );

	bool 	applyPionKinematicCuts( pion pi );
	bool 	applyPionKinematicCuts( genPion pi );
	bool 	applyPionDetectorCuts( pion pi, electron e );

	void	loadCutValues (int torusBending, double EBeam); //  -1 for In-bending, +1 for Out-bending
	bool	applyAcceptanceMatching( pion pi, int dim );

	void	printCutValues ();
	void	SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;}
	double	ComputeLightConeFraction ( TLorentzVector p );
	double	calcQStar ( TVector3 eP3, TVector3 piP3, double Ebeam );
	void loadMatchingFunctions( TString fileName );
	void loadMatchingFunctions();
	void loadAcceptanceMap( TString fileName );
	int acceptance_match_3d_cont( double phi_part, double theta, double p, int chargeIdx);
	int acceptance_match_3d( double phi_part, double theta, double p, int charge);
	int FindMatch(TVector3 p, clas12::mcpar_ptr mcparts, std::vector<int> part_list);
	
	
	bool checkElAcceptance( double p, double phi, double theta );
	bool checkPipAcceptance( double p, double phi, double theta );
	bool checkPimAcceptance( double p, double phi, double theta );
	
	//these are used in plotting program
	bool applyElectronFiducials( electron e );
	bool applyElectronPCAL( electron e );
	bool applyElectronEDep( electron e );
	bool applyElectronSF( electron e );
	bool applyElectronCorrelation( electron e );
	bool applyElectronVertex( electron e );
	bool applyPionDetectorFiducials( pion pi );
	bool applyPionDetectorChi2( pion pi );
	bool applyPionDetectorVertex( pion pi, electron e );
	
	bool checkAcceptance( double p, double phi, double theta, int particle );

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

	double acceptanceMap[6][10][3][3];  //sector, p bin, particle type, number of parameters	
	double fitBounds[6][5][3][3];	

};

#endif
