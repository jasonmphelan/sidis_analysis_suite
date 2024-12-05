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
	const int nRuns[3] = {7, 4, 3};
	const int monteCarloRuns[3][7] ={ 
				{7224, 7302, 7304, 7393, 7439, 7520, 8320},
				{ 8389, 7608, 7769, 8371, 8398, 0, 0 },
				{ 7427, 8376, 8436, 0, 0, 0, 0} 
				};

	//Rho subtraction of non exclusive background
	const double bac_norm = 0.178568;
	const double bac_min_mx = 1.15;
	const double bac_max_mx = 1.45; 

}

#endif
