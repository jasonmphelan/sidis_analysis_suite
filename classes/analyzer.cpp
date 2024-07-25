// Erez O. C., Oct-6, 2021
#include "analyzer.h"
#define CUT_PATH _DATA


analyzer::analyzer(int _fdebug_, int _torusBending_){
    //SetVerbosity    (_fdebug_);
    SetTorusBending (_torusBending_);
}

analyzer::~analyzer(){}

//bool 	applyElectronDetectorCuts( electron e );
//bool 	applyElectronKinematics( electron e );

//bool 	applyPionKinematics( pion pi );
//bool 	applyPionDetectorCuts( pion pi );

//void	loadAcceptanceMatching (int torusBending=1); //  -1 for In-bending, +1 for Out-bending
//void	loadCorrections (TFile corrFileName); //  -1 for In-bending, +1 for Out-bending
//void	getCorrections ( electron e, pion pi);
//void	applyAcceptanceMatching( pion pi );	
//void	printCutValues ();
//double	fillDetectorHistograms();



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



	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)

		int bending  = 1 ? (torusBending==-1) : 0;
		bool DC_fid  = dcfid.DC_fid_xy_sidis( 11,                 // particle PID,
						e_DC_x[regionIdx],  // x
						e_DC_y[regionIdx],  // y
						e.getDC_sector(),        // sector
						regionIdx+1,        // layer
						bending );           // torus bending
		if (DC_fid == false) { return false; }
	}	
	
	if( ! (// fiducial cuts on PCAL
		e.getW() > e_PCAL_W_min
		&&  e.getV() > e_PCAL_V_min)) return false;


	if( ! (
		// Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
		e.getEpcal() > e_E_PCAL_min)) return false;
		
	if(  !(epid.isElectron(&e)) ) return false;

	if( ! (
		// Sampling fraction cut
		( (e.getEpcal() + e.getEecin() + e.getEecout())/e.get3Momentum().Mag()) > SamplingFraction_min
		&& ( e.getEecin()/e.get3Momentum().Mag() > PCAL_ECIN_SF_min - e.getEpcal()/e.get3Momentum().Mag() ) // RGA AN puts "<" here mistakenly

		// Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
		// Spring 19 and Spring 2020 in-bending.
		// Fall	 2019 (without low-energy-run) was out-bending.
	) ) return false;
	
	/*
	if( !((aux.cutValue_Vz_min < v_e.Z()) && (v_e.Z() < aux.cutValue_Vz_max))
	) return false;
	*/

	if( ! ((e.getVt().Z() > -5) && (e.getVt().Z() < -1))
		) return false;
	
	
	//dcElectrons++;
	return true;
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

	
	int PDGcode;
	double    C;
	
	if (pi.getCharge() > 0){
		C       = 0.88;
	} 
	else if (pi.getCharge() < 0) {
		C       = 0.93;
	} 
	else {
		std::cout << "π charge ill-defined, returning false" << std::endl;
		return false;
	}

	double DC_x[3] = {pi.getDC_x1(), pi.getDC_x2(), pi.getDC_x3()};
	double DC_y[3] = {pi.getDC_y1(), pi.getDC_y2(), pi.getDC_y3()};
	double DC_z[3] = {pi.getDC_z1(), pi.getDC_z2(), pi.getDC_z3()};

	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
		
		int bending  = 1 ? (torusBending==-1) : 0;
		bool DC_fid  = dcfid.DC_fid_th_ph_sidis(pi.getPID(),            // particle PID
							DC_x[regionIdx],    // x
							DC_y[regionIdx],    // y
							DC_z[regionIdx],    // z
							pi.getDC_sector(),          // sector
							regionIdx+1,        // layer
							bending);           // torus bending
		
		if (DC_fid == false) { return false; }
	}


	//double chi2PID = pi.getDC_chi2() / pi.getDC_NDF();
	if(! (
	
	// pi+ Identification Refinement - chi2PID vs. momentum
	( Chi2PID_pion_lowerBound( pi.get3Momentum().Mag(), C ) < pi.getChi2()
         && pi.getChi2() < Chi2PID_pion_upperBound( pi.get3Momentum().Mag() , C ) )
	
	)) { return false; }
       
	//if( (cuts == 1 || cuts == 3 ) &&  !( fabs((v_e-v_pi).Z()) < aux.cutValue_Ve_Vpi_dz_max )
	//) { return false; }

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


bool analyzer::applyPionKinematicCuts( pion pi ){

	if ( ( pi.getMx() < Mx_min || pi.getMx() > Mx_max) ) { return false; }
	//if ( cut_type == 1 && ( M_x[i] < 1.5 || M_x[i] > 5.) ) { return false; }
	double p = pi.get3Momentum().Mag();
	if ( p < P_pi_min || p > P_pi_max ) { return false; }
	if ( pi.getZ() < Z_min  ||  pi.getZ() > Z_max ) { return false; }
	double theta = pi.get3Momentum().Theta();
	if ( theta < theta_min || theta > theta_max ){ return false; }
	return true;
}


bool analyzer::applyAcceptanceMatching( pion pi, int dim ){
	double theta = pi.get3Momentum().Theta()*rad_to_deg;
	double p = pi.get3Momentum().Mag();
	int sector_i = pi.getDC_sector();
	
	if( dim == 2 ){
		
		double acc_map_pip_min = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/p;                      
                double acc_map_pim_min = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/p;

		if ( theta > acc_map_pip_min && theta > acc_map_pim_min ){return true;}
		else { return false; }
	}
	//else if( dim == 3 ){

	//	double phi = pi.get3Momentum().Phi()*rad_to_deg;
		//if( acceptance_match_3d( phi, theta, p, 0 ) && acceptance_match_3d( phi, theta, p, 1) ){
	//	if( acceptance_match_3d_cont( phi, theta, p, match3d ) ){ return true; }
	//	else{ return false; }
	
	//}
	else{
		std::cout<<"Bad argument for dimensionality of acceptance matching... returning false\n";
		return false;
	}
}

bool acceptance_match_3d( double phi_part, double theta, double p, int charge){
	
	//set momentum bins
	int this_bin_p;
	for( int i = 0; i < 4; i++ ){
		if( p > p_bin_edges[i] && p < p_bin_edges[i+1]){ this_bin_p = i; }
	}

	bool passCut = false;
	
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

		if( theta > theta_min ){ passCut = true; }		
	}
	return passCut;
}
/*
bool acceptance_match_3d_cont( double phi_part, double theta, double p, TF1 * fitFuncs[6][3]){
	
	//set momentum bins
	//int this_bin_p;
	//for( int i = 0; i < 4; i++ ){
	//	if( p > p_bin_edges[i] && p < p_bin_edges[i+1]){ this_bin_p = i; }
	//}

	bool passCut[2] = {false, false};
	for( int charge = 1; charge <= 2; charge++){
		for( int sector = 0; sector < 6; sector++ ){
			double phi = phi_part;
			if( sector ==3 && phi < 100. ){ phi += 360; }

			//Get parameters from constants
			double theta_0 = fitFuncs[sector][0]->Eval(p);//phi_theta_bowl_theta_min[sector - 1][this_bin_p][1];
			double phi_0 = fitFuncs[sector][charge]->Eval(p);//phi_theta_bowl_phi0[sector-1][this_bin_p][charge];

			//compute cut value
			double theta_min = theta_0 + pow( (phi-phi_0), 2 )/( theta_bowl_width - pow( ( phi - phi_0 ), 2 ) );

			//If phi is outside bowl, set theta_min = theta_max
			if( theta_bowl_width - pow( (phi - phi_0), 2 ) < 0 ){
				theta_min = 35.;
			}
			//cout<<"theta min cont"<<theta_min<<endl;
			if( theta > theta_min ){ passCut[charge-1] = true; }		
		}
	}
	return (passCut[0] && passCut[1]);
}
*/
/*
void analyzer::loadCutValues(int torusBending){
// read cut values csv file

char cutFileName[100];
sprintf(cutFileName,"%s/BANDcutValues.csv",std::string(_DATA).c_str());

cutValues = csvr.read_csv(cutFileName);

if (torusBending==-1){ // in-bending torus field
// Spring 19 and Spring 2020 in-bending.
cutValue_Vz_min = FindCutValue("Vz_e_min_inbending");
cutValue_Vz_max = FindCutValue("Vz_e_max_inbending");
} else if (torusBending==1){ // Out-bending torus field
// Fall 2019 (without low-energy-run) was out-bending.
cutValue_Vz_min = FindCutValue("Vz_e_min_outbending");
cutValue_Vz_max = FindCutValue("Vz_e_max_outbending");

    } else {
        std::cout
        << "Un-identified torus bending "
        << torusBending
        << ", return" << std::endl;
        return;
    }
    
    cutValue_e_PCAL_W               = FindCutValue("e_PCAL_W_min");
    cutValue_e_PCAL_V               = FindCutValue("e_PCAL_V_min");
    cutValue_e_E_PCAL               = FindCutValue("e_E_PCAL_min");
    cutValue_SamplingFraction_min   = FindCutValue("SamplingFraction_min");
    cutValue_PCAL_ECIN_SF_min       = FindCutValue("PCAL_ECIN_SF_min");
    cutValue_Ve_Vpi_dz_max          = FindCutValue("(Ve-Vpi)_z_max");
    cutValue_Q2_min                 = FindCutValue("Q2_min");
    cutValue_Q2_max                 = FindCutValue("Q2_max");
    cutValue_W_min                  = FindCutValue("W_min");
    cutValue_y_max                  = FindCutValue("y_max");
    cutValue_e_theta_min            = FindCutValue("e_theta_min");
    cutValue_e_theta_max            = FindCutValue("e_theta_max");
    cutValue_pi_theta_min           = FindCutValue("pi_theta_min");
    cutValue_pi_theta_max           = FindCutValue("pi_theta_max");
    cutValue_Ppi_min                = FindCutValue("Ppi_min");
    cutValue_Ppi_max                = FindCutValue("Ppi_max");
    cutValue_Pe_min                 = FindCutValue("Pe_min");
    cutValue_Pe_max                 = FindCutValue("Pe_max");
    cutValue_Zpi_min                = FindCutValue("Zpi_min");
    cutValue_Zpi_max                = FindCutValue("Zpi_max");
    
    if (fdebug>2) { printCutValues(); }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void analyzer::printCutValues(){
    std::cout << "Using the following cut values:" << std::endl;
    std::cout <<
    "Vz_min: "                 << cutValue_Vz_min                 << ", " << std::endl <<
    "Vz_max: "                 << cutValue_Vz_max                 << ", " << std::endl <<
    "e_PCAL_W: "               << cutValue_e_PCAL_W               << ", " << std::endl <<
    "e_PCAL_V: "               << cutValue_e_PCAL_V               << ", " << std::endl <<
    "e_E_PCAL: "               << cutValue_e_E_PCAL               << ", " << std::endl <<
    "SamplingFraction_min: "   << cutValue_SamplingFraction_min   << ", " << std::endl <<
    "PCAL_ECIN_SF_min: "       << cutValue_PCAL_ECIN_SF_min       << ", " << std::endl <<
    "Ve_Vpi_dz_max: "          << cutValue_Ve_Vpi_dz_max          << ", " << std::endl <<
    "Q2_min: "                 << cutValue_Q2_min                 << ", " << std::endl <<
    "Q2_max: "                 << cutValue_Q2_max                 << ", " << std::endl <<
    "W_min: "                  << cutValue_W_min                  << ", " << std::endl <<
    "y_max: "                  << cutValue_y_max                  << ", " << std::endl <<
    "e_theta_min: "            << cutValue_e_theta_min            << ", " << std::endl <<
    "e_theta_max: "            << cutValue_e_theta_max            << ", " << std::endl <<
    "pi_theta_min: "           << cutValue_pi_theta_min           << ", " << std::endl <<
    "pi_theta_max: "           << cutValue_pi_theta_max           << ", " << std::endl <<
    "Ppi_min: "                << cutValue_Ppi_min                << ", " << std::endl <<
    "Ppi_max: "                << cutValue_Ppi_max                << ", " << std::endl <<
    "Pe_min: "                 << cutValue_Pe_min                 << ", " << std::endl <<
    "Pe_max: "                 << cutValue_Pe_max                 << ", " << std::endl <<
    "Zpi_min: "                << cutValue_Zpi_min                << ", " << std::endl <<
    "Zpi_max: "                << cutValue_Zpi_max               << ", " << std::endl <<
    std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void analyzer::SetTorusBendingFromRunNumber ( Int_t RunNumber ){
    // -1 for In-bending, +1 for Out-bending
    // For BAND data
    // Spring 19 and Spring 2020 was in-bending
    // Fall 2019 (without low-energy-run) was out-bending
    
    if (6420 <= RunNumber && RunNumber <= 6598){
        this->torusBending = -1;
    }
    else if (11362 <= RunNumber && RunNumber <= 11571){
        this->torusBending = -1;
    }
    else if (6164 <= RunNumber && RunNumber <= 6399){
        this->torusBending = +1;
    }
    else{
        this->torusBending = 0;
    }
}
*/


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double analyzer::ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    double m = p.Mag();
    double alpha = (p.E() - p.Z())/m;
    return alpha;
}

/*
// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool analyzer::eepiPassedKinematicalCriteria(Double_t Ebeam,
                                                          Double_t omega,
                                                          Double_t Q2,
                                                          Double_t y,
                                                          Double_t W,
                                                          TLorentzVector pi,
                                                          TLorentzVector e){
    
    if(   (      cutValue_Q2_min < Q2             &&             Q2 < cutValue_Q2_max       )
       && (       cutValue_W_min < W                                                        )
       && (                                                       y < cutValue_y_max        )
       && ( cutValue_e_theta_min < e.Theta()*r2d  &&  e.Theta()*r2d < cutValue_e_theta_max  )
       && (cutValue_pi_theta_min < pi.Theta()*r2d && pi.Theta()*r2d < cutValue_pi_theta_max )
       && (     cutValue_Ppi_min < pi.P()         &&         pi.P() < cutValue_Ppi_max      )
       && (      cutValue_Pe_min < e.P()          &&          e.P() < Ebeam                 )
       && (     cutValue_Zpi_min < Zpi            &&            Zpi < cutValue_Zpi_max      )
       ) {
        if (fdebug>2) {
            std::cout << "succesfully passed (e,e'π) kinematical cuts"
            << std::endl
            << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;
        }
        return true;
    }
    
    return false;
}
*/


// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double analyzer::calcQStar(TVector3 eP3, TVector3 piP3, double Ebeam){
    //
    // Compute the novel SIDIS observable, q∗
    // designed to be maximally resilient against resolution
    // effects while delivering the same sensitivity to TMD dynamics as Pt
    //
    // [Gao et al., "A Better Angle on Hadron Transverse Momentum Distributions at the EIC"]
    // arXiv:2209.11211v1
    //
    // by Natalie Wright Oct-12, 2022
    //
    // input
    // ------
    // eP3      TVector3    electron momentum (in q-Frame)
    // piP3     TVector3    π momentum (in q-Frame)
    //
    // return
    // ------
    // qstar    double      coplanarity momentum transfer part
    //
    
    
    
    double tan_phi_acop = piP3.Y()/piP3.X(); // EIC Frame
    double       eta_pi = TMath::ATanH(-piP3.Z()/piP3.Mag()); // - sign from EIC frame
    double        eta_e = TMath::ATanH(-eP3.Z()/eP3.Mag());
    double    delta_eta = eta_pi - eta_e;
    
    double        qstar =   2 * Ebeam * TMath::Exp(eta_pi) * tan_phi_acop
                            /
                            (1 + TMath::Exp(delta_eta));
    
    return qstar;    
}



