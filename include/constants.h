#ifndef CONSTANTS_H_
#define CONSTANTS_H_
#include "TLorentzVector.h"

namespace constants{
	//Useful Constants
	const double rad_to_deg = 180./TMath::Pi();

	//Particle Infor
	const double	Me  = 0.00051099895; // GeV/c2
	const double	Mpi = 0.139570; // GeV/c2
	const double	Mp  = 0.938272; // GeV/c2
	const double	Mn = 0.939565; // GeV/c2
	const double	Md = 1.875; // GeV/c2
	
	const TLorentzVector p_rest( 0, 0, 0, Mp );
	const TLorentzVector d_rest( 0, 0, 0, Md );

	//Detector info
	const int DC_layers[3] = {6, 18, 36};
	const double n_rich = 1.05;
	
	//Monte Carlo info
	const int nRuns[3] = {6, 3, 1};
	const int monteCarloRuns[3][6] ={ 
				{7224, 7302, 7304, 7393, 7439, 7520},
				{ 7607, 7608, 7769 },
				{ 7427 } 
				};
}

#endif
