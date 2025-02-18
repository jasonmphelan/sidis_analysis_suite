#include "analyzer.h"
#include "TFile.h"
#include "TF1.h"
#define CUT_PATH _DATA


analyzer::analyzer(int _fdebug_, int _torusBending_){
    SetTorusBending (_torusBending_);
}

analyzer::~analyzer(){}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double analyzer::Chi2PID_pion_lowerBound( Double_t p, Double_t C){
    // compute lower bound for chi2PID for a pi+
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 75
    // "Strict cut"
    //
    // input:
    // -------
    // p        pion momentum
    // C        is the scaling factor for sigma (away from the mean value of the pion distribution)
    //
    
    return ( -C * 3 );
}

// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double analyzer::Chi2PID_pion_upperBound( Double_t p, Double_t C){
    // compute upper bound for chi2PID for a pi+
    // based on RGA_Analysis_Overview_and_Procedures_Nov_4_2020-6245173-2020-12-09-v3.pdf
    // p. 75
    // "Strict cut"
    //
    // input:
    // -------
    // p        pion momentum
    // C        is the scaling factor for sigma (away from the mean value of the pion distribution)
    //
    
    if (p<2.44)
        
        return C*3;
    
    else if (p<4.6)
        
        return C*( 0.00869 + 14.98587*exp(-p/1.18236)+1.81751*exp(-p/4.86394) ) ;
    
    else
        
        return C*( -1.14099 + 24.14992*exp(-p/1.36554) + 2.66876*exp(-p/6.80552) );
    
}

void analyzer::loadCutValues(int torusBending, double EBeam){
	epid.setParamsRGB(EBeam);
}


bool analyzer::applyElectronDetectorCuts( electron e ){
	if (e.getDC_sector() == 0) return false;


	double e_DC_x[3] = {e.getDC_x1(), e.getDC_x2(), e.getDC_x3()};
	double e_DC_y[3] = {e.getDC_y1(), e.getDC_y2(), e.getDC_y3()};
	double e_DC_z[3] = {e.getDC_z1(), e.getDC_z2(), e.getDC_z3()};


/*
	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
		if( e.getEdge(regionIdx) < e_fid_cuts[regionIdx] ){ return false; }
		
		//int bending  = 1 ? (torusBending==-1) : 0;
		//bool DC_fid  = dcfid.DC_fid_xy_sidis( 11,                 // particle PID,
		//				e_DC_x[regionIdx],  // x
		//				e_DC_y[regionIdx],  // y
		//				e.getDC_sector(),        // sector
		//				regionIdx+1,        // layer
		//				bending );           // torus bending
		//if (DC_fid == false) { return false; }
		
	}	
*/
	// PCAL FIDUCIAL	
	if( ! (	e.getW() > e_PCAL_W_min	&&  e.getV() > e_PCAL_V_min)) return false;


	//PCAL MIN EDEP CUT
	if( ! ( e.getEpcal() > e_E_PCAL_min) ) return false;
		
	//Electron SF cut
	if(  !(epid.isElectron(&e)) ) return false;

	//if( ! (( (e.getEpcal() + e.getEecin() + e.getEecout())/e.get3Momentum().Mag()) > SamplingFraction_min ) ) return false;
	
	//SF CORRELATION CUT
	if( !( e.getEecin()/e.get3Momentum().Mag() > PCAL_ECIN_SF_min - e.getEpcal()/e.get3Momentum().Mag() )) return false;
	
	//ELECTRON VERTEX CUT
	if( ! ((e.getVt().Z() > Vz_e_min_inbending) && (e.getVt().Z() < Vz_e_max_inbending))) return false;
	
	
	return true;
}
bool analyzer::applyElectronFiducials( electron e ){
	if (e.getDC_sector() == 0) return false;


	double e_DC_x[3] = {e.getDC_x1(), e.getDC_x2(), e.getDC_x3()};
	double e_DC_y[3] = {e.getDC_y1(), e.getDC_y2(), e.getDC_y3()};
	double e_DC_z[3] = {e.getDC_z1(), e.getDC_z2(), e.getDC_z3()};

	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
		if( e.getEdge(regionIdx) < e_fid_cuts[regionIdx] ){ return false; }
	}	
	return true;
}

bool analyzer::applyElectronPCAL( electron e ){
	return ( e.getW() > e_PCAL_W_min	&&  e.getV() > e_PCAL_V_min);
}

bool analyzer::applyElectronEDep( electron e ){

	//PCAL MIN EDEP CUT
	return  ( e.getEpcal() > e_E_PCAL_min) ;
}
bool analyzer::applyElectronSF( electron e ){
	//Electron SF cut
	return (epid.isElectron(&e)) ;
}
	
bool analyzer::applyElectronCorrelation( electron e ){
	//SF CORRELATION CUT
	return(  e.getEecin()/e.get3Momentum().Mag() > PCAL_ECIN_SF_min - e.getEpcal()/e.get3Momentum().Mag() );
}


bool analyzer::applyElectronVertex( electron e ){
	//ELECTRON VERTEX CUT
	return ((e.getVt().Z() > -5) && (e.getVt().Z() < -1));
}

bool analyzer::applyPionDetectorCuts( pion pi, electron e ){
	// decide if pion (pi+ or pi-) passed event selection cuts
	//
	// input:
	// --------
	// DC_x, DC_y   pi drift-chamber coordinates
	// chi2PID      pi chi2PID     (pips_chi2PID)
	// p            pi momentum    (pi.P())
	//

	if (pi.getDC_sector() == 0) { return false;}

	double    C;
	
	if (pi.getCharge() > 0){ C = 0.88;} 
	else if (pi.getCharge() < 0) {C = 0.93; } 
	else {
		std::cout << "π charge ill-defined, returning false" << std::endl;
		return false;
	}

	//double DC_x[3] = {pi.getDC_x1(), pi.getDC_x2(), pi.getDC_x3()};
	//double DC_y[3] = {pi.getDC_y1(), pi.getDC_y2(), pi.getDC_y3()};
	//double DC_z[3] = {pi.getDC_z1(), pi.getDC_z2(), pi.getDC_z3()};
/*
	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
	
		if( pi.getEdge(regionIdx) < pi_fid_cuts[ (int)(pi.getCharge() < 0 ) ][regionIdx] ){ return false; }
		
		//int bending  = 1 ? (torusBending==-1) : 0;
		//bool DC_fid  = dcfid.DC_fid_th_ph_sidis(pi.getPID(),            // particle PID
		//					DC_x[regionIdx],    // x
		//					DC_y[regionIdx],    // y
		//					DC_z[regionIdx],    // z
		//					pi.getDC_sector(),          // sector
		//					regionIdx+1,        // layer
		//					bending);           // torus bending
		//
		//if (DC_fid == false) { return false; }
	
	}
	
*/
	

	//PION CHI2 vs P CUT
	if(! (	( Chi2PID_pion_lowerBound( pi.get3Momentum().Mag(), C ) < pi.getChi2()
         	&& pi.getChi2() < Chi2PID_pion_upperBound( pi.get3Momentum().Mag() , C ) ))) 
		{return false; }
       
	//DELTA VERTEX CUT
	if( !( (pi.getVt() - e.getVt()).Z() > -7 && (pi.getVt() - e.getVt()).Z() < 5 ) ) { return false; }
	//if( ! ( abs( (pi.getVt() - e.getVt()).Z() ) < 20 ) ) { return false; }
	return true;
}

bool analyzer::applyPionDetectorFiducials( pion pi ){
	// decide if pion (pi+ or pi-) passed event selection cuts
	//
	// input:
	// --------
	// DC_x, DC_y   pi drift-chamber coordinates
	// chi2PID      pi chi2PID     (pips_chi2PID)
	// p            pi momentum    (pi.P())
	//

	

	double DC_x[3] = {pi.getDC_x1(), pi.getDC_x2(), pi.getDC_x3()};
	double DC_y[3] = {pi.getDC_y1(), pi.getDC_y2(), pi.getDC_y3()};
	double DC_z[3] = {pi.getDC_z1(), pi.getDC_z2(), pi.getDC_z3()};

	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
		if( pi.getEdge(regionIdx) < pi_fid_cuts[ (int)(pi.getCharge() < 0 ) ][regionIdx] ){ 
			return false; 
		}
	}
	return true;
}
bool analyzer::applyPionDetectorChi2( pion pi ){
	double    C;
	
	if (pi.getCharge() > 0){ C = 0.88;} 
	else if (pi.getCharge() < 0) {C = 0.93; } 
	else {
		std::cout << "π charge ill-defined, returning false" << std::endl;
		return false;
	}
	//PION CHI2 vs P CUT
	return (	( Chi2PID_pion_lowerBound( pi.get3Momentum().Mag(), C ) < pi.getChi2()
         	&& pi.getChi2() < Chi2PID_pion_upperBound( pi.get3Momentum().Mag() , C ) )); 
}

bool analyzer::applyPionDetectorVertex( pion pi, electron e ){
       
	//DELTA VERTEX CUT
	if( !( (pi.getVt() - e.getVt()).Z() > -7 && (pi.getVt() - e.getVt()).Z() < 5 ) ) { return false; }
	
	return true;
}



bool analyzer::applyElectronKinematicCuts( electron e ){
		if( sqrt(e.getW2()) < W_min ) { return false; }
		if( e.getQ2() < Q2_min || e.getQ2() > Q2_max ) { return false; }
                if( e.getXb() < xB_min || e.getXb() > xB_max ) { return false; }
                if( e.getY() > y_max ) { return false; }
                
		double p = e.get3Momentum().Mag();
		if( p < P_e_min || p > P_e_max ) { return false; }
                
		double theta = e.get3Momentum().Theta()*rad_to_deg;
		if( theta < theta_min || theta > theta_max ){ return false; }
		
		return true;
}
bool analyzer::applyElectronKinematicCuts( genElectron e ){
		if( sqrt(e.getW2()) < W_min ) { return false; }
		if( e.getQ2() < Q2_min || e.getQ2() > Q2_max ) { return false; }
                if( e.getXb() < xB_min || e.getXb() > xB_max ) { return false; }
                if( e.getY() > y_max ) { return false; }
                
		double p = e.get3Momentum().Mag();
		if( p < P_e_min || p > P_e_max ) { return false; }
                
		double theta = e.get3Momentum().Theta()*rad_to_deg;
		if( theta < theta_min || theta > theta_max ){ return false; }
		
		return true;
}


bool analyzer::applyPionKinematicCuts( pion pi ){


	if ( ( pi.getMx() < Mx_min || pi.getMx() > Mx_max) ) { return false; }
	double p = pi.get3Momentum().Mag();
	if ( p < P_pi_min || p > P_pi_max ) { return false; }
	if ( pi.getZ() < Z_min  ||  pi.getZ() > Z_max ) { return false; }
	double theta = pi.get3Momentum().Theta()*rad_to_deg;
	if ( theta < theta_min || theta > theta_max ){ return false; }
	return true;
}
bool analyzer::applyPionKinematicCuts( genPion pi ){


	if ( ( pi.getMx() < Mx_min || pi.getMx() > Mx_max) ) { return false; }
	double p = pi.get3Momentum().Mag();
	if ( p < P_pi_min || p > P_pi_max ) { return false; }
	if ( pi.getZ() < Z_min  ||  pi.getZ() > Z_max ) { return false; }
	double theta = pi.get3Momentum().Theta()*rad_to_deg;
	if ( theta < theta_min || theta > theta_max ){ return false; }
	return true;
}

//Function to do all acceptance matching
bool analyzer::applyAcceptanceMatching( pion pi, int dim ){
	double theta = pi.get3Momentum().Theta()*rad_to_deg;
	double p = pi.get3Momentum().Mag();
	int sector_i = pi.getDC_sector();
	
	if( dim == 2 ){
		return acceptance_match_2d( theta, p, sector_i );	
	}
	else if( dim == 3 ){

		double phi = pi.get3Momentum().Phi()*rad_to_deg;
		//if( acceptance_match_3d( phi, theta, p, 0 ) && acceptance_match_3d( phi, theta, p, 1) ){
		if( acceptance_match_3d_cont( phi, theta, p, 0 ) > -1
			&& acceptance_match_3d_cont( phi, theta, p, 1 ) > -1   ){ return true; }
		else{ return false; }
	
	}
	else{
		std::cout<<"Bad argument for dimensionality of acceptance matching... returning false\n";
		return false;
	}
}

//2d acceptance matching cut
bool analyzer::acceptance_match_2d( double theta, double p, int sector_i){
		//double acc_map_pip_min = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/p;                      
                //double acc_map_pim_min = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/p;

		bool max_pip = (bool) ( theta < match2d[sector_i - 1][0][0]->Eval(p) );
		bool max_pim = (bool) ( theta < match2d[sector_i - 1][0][1]->Eval(p) );
		
		bool min_pip = (bool) ( theta > match2d[sector_i - 1][1][0]->Eval(p) );
		bool min_pim = (bool) ( theta > match2d[sector_i - 1][1][1]->Eval(p) );
		
		return (max_pip && max_pim && min_pip && min_pim);

		//return (theta > acc_map_pip_min && theta > acc_map_pim_min );
}

//Discrete 3d accepance matching cut
int analyzer::acceptance_match_3d( double phi_part, double theta, double p, int charge){
	
	//set momentum bins
	int this_bin_p;
	for( int i = 0; i < 4; i++ ){
		if( p > p_bin_edges[i] && p < p_bin_edges[i+1]){ this_bin_p = i; }
	}

	int passCut = -1;
	
	for( int sector = 1; sector <=6; sector++ ){
		double phi = phi_part;
		if( sector ==4 && phi < 100. ){ phi += 360; }

		//Get parameters from constants
		double theta_0 = phi_theta_bowl_theta_min[sector - 1][this_bin_p][1];
		double phi_0 = phi_theta_bowl_phi0[sector-1][this_bin_p][charge];

		//compute cut value
		double theta_min = theta_0 + pow( (phi-phi_0), 2 )/( theta_bowl_width - pow( ( phi - phi_0 ), 2 ) );

		//If phi is outside bowl, set theta_min = theta_max
		if( theta_bowl_width - pow( (phi - phi_0), 2 ) < 0 ){
			theta_min = 35.;
		}

		if( theta > theta_min ){ passCut = sector - 1; }		
	}
	return passCut;
}

//Continuous acceptance matching cut
int analyzer::acceptance_match_3d_cont( double phi_part, double theta, double p, int chargeIdx){
	if( chargeIdx < 0 || chargeIdx > 1 ){
		std::cout<<"Invalid Charge Index... returning -1\n";
		return -1;
	}
	
	for( int sector = 0; sector < 6; sector++ ){
		double phi = phi_part;
		if( sector ==3 && phi < 100. ){ phi += 360; }
		
		//Get parameters from functions
		double theta_0 = match3d[sector][0]->Eval(p);//phi_theta_bowl_theta_min[sector - 1][this_bin_p][1];
		double phi_0 = match3d[sector][chargeIdx+1]->Eval(p);//phi_theta_bowl_phi0[sector-1][this_bin_p][charge];

		//compute cut value
		double theta_min = theta_0 + pow( (phi-phi_0), 2 )/( theta_bowl_width - pow( ( phi - phi_0 ), 2 ) );

		//If phi is outside bowl, set theta_min = theta_max
		if( theta_bowl_width - pow( (phi - phi_0), 2 ) < 0 ){
			theta_min = 35.;
		}
		if( theta > theta_min ){ return sector; }		
	}
	return -1;
}


//Load matching functions
void analyzer::loadMatchingFunctions( TString fileName ){

	TFile matchFile2D( (TString) CUT_PATH + "/acceptance_matching/" + fileName);

	TString chargeType[2] = {"pip", "pim"};	
	TString boundType[2] = {"max", "min"};
	for( int i = 0; i < 6; i++ ){		
		for( int j = 0; j < 2; j++ ){
			for( int k = 0; k < 2; k++ ){
				match2d[i][j][k] = (TF1 *)matchFile2D.Get(boundType[j] + Form("_%i_", i) + chargeType[k]);
			}
		}
	}
}

void analyzer::loadMatchingFunctions3D( TString fileName ){

	TFile matchFile3D( (TString) CUT_PATH + "/acceptance_matching/" + fileName);

	TString chargeType[2] = {"pip", "pim"};	

	for( int i = 0; i < 6; i++ ){		
		match3d[i][0] = (TF1 *)matchFile3D.Get(Form("fTheta0_%i", i));
		match3d[i][1] = (TF1 *)matchFile3D.Get(Form("fPhi0_pip_%i", i));
		match3d[i][2] = (TF1 *)matchFile3D.Get(Form("fPhi0_pim_%i", i));
	}
}

void analyzer::loadMatchingFunctions(){
	loadMatchingFunctions("matchCut2D.root");
}
void analyzer::loadMatchingFunctions3D(){
	loadMatchingFunctions("matchCut3D.root");
}


//Find gemc/gen match
int analyzer::FindMatch(TVector3 p, clas12::mcpar_ptr mcparts, std::vector<int> part_list){

	double temp_min_dp = 9999;
	double p_idx = -1;
	double temp_min_dtheta = 9999;
	double theta_idx = -1;
	double temp_min_dphi = 9999;
	double phi_idx = -1;

	for ( int i : part_list ){
		mcparts->setEntry(i);

		if( abs( p.Mag() - mcparts->getP() ) < temp_min_dp ){	
			temp_min_dp = abs( p.Mag() - mcparts->getP() );
			p_idx = i;
		}		
		if( abs( p.Theta() - mcparts->getTheta() ) < temp_min_dtheta ){	
			temp_min_dtheta = abs( p.Theta() - mcparts->getTheta() );
			theta_idx = i;
		}	
		if( abs( p.Phi() - mcparts->getPhi() ) < temp_min_dphi ){	
			temp_min_dphi = abs( p.Phi() - mcparts->getPhi() );
			phi_idx = i;
		}	
	}
	
	if( theta_idx != p_idx || p_idx != phi_idx){
		return -1;
	}
	
	return p_idx;
}

//Load acceptance map
void analyzer::loadAcceptanceMap(TString fileName){
	TFile f(fileName);

	TString parType[3] = {"e", "pip", "pim"};
	int nPbins[3] = {10, 5, 5};
	for( int par = 0; par < 3; par++ ){
		for( int p = 0; p < nPbins[par]; p++ ){
			for( int sec = 0; sec < 6; sec++ ){
				TF1 * f1 = (TF1 *) f.Get(Form("fThetaPhi_sec_%i_bin_%i_%s", sec, p, parType[par].Data() ) );
				TVector3 * fBounds = (TVector3 *) f.Get(Form("fitBounds_sec_%i_bin_%i_%s", sec, p, parType[par].Data() ) );
				
				acceptanceMap[sec][p][par][0] = f1->GetParameter(0);
				acceptanceMap[sec][p][par][1] = f1->GetParameter(1);
				acceptanceMap[sec][p][par][2] = f1->GetParameter(2);
				
				fitBounds[sec][p][par][0] = fBounds->X();
				fitBounds[sec][p][par][1] = fBounds->Y();
				fitBounds[sec][p][par][2] = fBounds->Z();
			}
		}
	}
}

//Check if in acceptance map
int analyzer::checkAcceptance( double p, double phi, double theta, int particle ){
	double par_0 = -999;
	double par_1 = -999;
	double par_2 = -999;

	double max = -999;
	double lower = 999;
	double upper = -999;

	//acceptanceMap[6][10][3][3];  //sector, p bin, particle type, number of parameters	
	double p_max = 5 + 5*( (int) (particle < 1) );
	int p_bin = (int)( ( (p)/(p_max) )*p_max);


	double cutMin = 40;
	for( int sec = 0; sec < 6; sec++ ){
		double phi_temp = phi;
		par_0 = acceptanceMap[sec][p_bin][particle][0];
		par_1 = acceptanceMap[sec][p_bin][particle][1];
		par_2 = acceptanceMap[sec][p_bin][particle][2];
		
		max = fitBounds[sec][p_bin][particle][2];
		lower = fitBounds[sec][p_bin][particle][0];
		upper = fitBounds[sec][p_bin][particle][1];
	
		if( sec ==3 && phi_temp < 100. ){ phi_temp += 360; }
		cutMin = par_0*(phi_temp - par_1)*(phi_temp - par_1) + par_2;
		
		if( theta > cutMin && theta < 40
			&& theta < max && phi_temp > lower && phi_temp < upper ){ 
			
			return sec; 
		}
	}

	return -1;
}

