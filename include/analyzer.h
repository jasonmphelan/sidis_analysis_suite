#ifndef __ANALYZER__
#define __ANALYZER__

#include "csv_reader.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
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

	void printCuts();

	void loadSamplingFractionParams( TString sfFile_name );
	void loadSamplingFractionParams( );

	void randomizeCuts();
	void	loadCutValues (int torusBending, double EBeam); //  -1 for In-bending, +1 for Out-bending
	bool	applyAcceptanceMatching( pion pi, int dim );

	void	printCutValues ();
	void	SetTorusBending (int _torusBending_)   {torusBending = _torusBending_;}
	double	ComputeLightConeFraction ( TLorentzVector p );
	double	calcQStar ( TVector3 eP3, TVector3 piP3, double Ebeam );
	void loadMatchingFunctions( TString fileName );
	void loadMatchingFunctions();
	void loadMatchingFunctions3D( TString fileName );
	void loadMatchingFunctions3D();
	void loadAcceptanceMap( TString fileName );
	int acceptance_match_3d_cont( double phi_part, double theta, double p, int chargeIdx);
	int acceptance_match_3d( double phi_part, double theta, double p, int charge);
	bool acceptance_match_2d( double theta, double p, int sector_i);
	int FindMatch(TVector3 p, clas12::mcpar_ptr mcparts, std::vector<int> part_list);
	
	
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
	
	int checkAcceptance( double p, double phi, double theta, int particle );
	
private:
	int	fdebug;
	int	torusBending; // -1 for In-bending, +1 for Out-bending
	int 	mode;


	TH3F *	accCorrection;
	TH3F *	binCorrection;
	TH3F *	kaonCorrection;
	TH3F *	rhoCorrection;

	TF1 * match3d[6][3];
	TF1 * match2d[6][2][2]; //[sector][max/min][pip/pim]

	e_pid epid;
	DCfid_SIDIS dcfid;
	electron * e;
	pion * pi;

	double SF_p_mean[6][3];
	double SF_p_sigma[6][3];

	double acceptanceMap[6][10][3][3];  //sector, p bin, particle type, number of parameters	
	double fitBounds[6][10][3][3];	

	//parameter modifiers for cut sensitivity
	double mod_el_fid[3] = {1, 1, 1};
	double mod_el_PCAL[2] = {1, 1};
	double mod_el_Edep = 1;
	double mod_SF_sigma = 1;
	double mod_el_corr = 1;
	double mod_el_vtz[2] = {1, 1};

	double mod_pi_fid[2][3] = {0};
	double mod_pi_vtz =  1;

	double mod_W = 1;
	double mod_Q2[2] = {1,1};
	double mod_xB[2] = {1, 1};
	double mod_y = 1;
	double mod_pe[2] = {1, 1};

	double mod_Mx = 1;
	double mod_ppi[2] = {1, 1};
	std::random_device rd;
};

#endif
