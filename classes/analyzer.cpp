#include "analyzer.h"
#include "TFile.h"
#include "TF1.h"
#include <random>
#define CUT_PATH _DATA


analyzer::analyzer(int _fdebug_, int _torusBending_){
	setAnalyzerLevel(mode);
	SetTorusBending (_torusBending_);
	std::fill_n( mod_pi_fid[0], 3, 1 );
	std::fill_n( mod_pi_fid[1], 3, 1 );
	//loadSamplingFractionParams();
}

analyzer::~analyzer(){}

void analyzer::printCuts(){
	std::cout<<"**************** Electron Selection Cuts ******************\n";
	for( int i = 0; i < 3; i++ ){
		std::cout<<"Electron Fiducial Region "<<i<<" : "<<mod_el_fid[i]*e_fid_cuts[i]<<std::endl;
	}
	std::cout<<"PCAL W : "<<mod_el_PCAL[0]*e_PCAL_W_min<<std::endl;
	std::cout<<"PCAL V : "<<mod_el_PCAL[1]*e_PCAL_V_min<<std::endl;
	std::cout<<"Electron Min Vertex : "<<mod_el_vtz[0]*Vz_e_min_inbending<<std::endl;
	std::cout<<"Electron Max Vertex : "<<mod_el_vtz[1]*Vz_e_max_inbending<<std::endl;
	std::cout<<"Minimum Edep : "<<mod_el_Edep*e_E_PCAL_min<<std::endl;
	std::cout<<"SF Sigma : "<<mod_SF_sigma*3.5<<std::endl;	
	std::cout<<"SF Correlation : "<<mod_el_corr*PCAL_ECIN_SF_min<<std::endl;

	std::cout<<"**************** Pion Selection Cuts ******************\n";
	for( int i = 0; i < 3; i++ ){
		std::cout<<"Pip Fiducial Region "<<i<<" : "<<mod_pi_fid[0][i]*pi_fid_cuts[0][i]<<std::endl;
		std::cout<<"Pim Fiducial Region "<<i<<" : "<<mod_pi_fid[1][i]*pi_fid_cuts[1][i]<<std::endl;
	}
	std::cout<<"Pip-e Min Vertex : "<<Vz_pi_mean[mode][0] - 3.5*mod_pi_vtz*Vz_pi_sigma[mode][0]<<std::endl;
	std::cout<<"Pip-e Max Vertex : "<<Vz_pi_mean[mode][0] + 3.5*mod_pi_vtz*Vz_pi_sigma[mode][0]<<std::endl;

	std::cout<<"Pim-e Min Vertex : "<<Vz_pi_mean[mode][1] - 3.5*mod_pi_vtz*Vz_pi_sigma[mode][1]<<std::endl;
	std::cout<<"Pim-e Max Vertex : "<<Vz_pi_mean[mode][1] + 3.5*mod_pi_vtz*Vz_pi_sigma[mode][1]<<std::endl;

	
	
	std::cout<<"**************** Kinematical Cuts ******************\n";
	std::cout<<"Minimum pe : "<<mod_pe[0]*P_e_min<<std::endl;
	std::cout<<"Maximum pe : "<<mod_pe[1]*P_e_max<<std::endl;
	std::cout<<"Minimum p_pi : "<<mod_ppi[0]*P_pi_min<<std::endl;
	std::cout<<"Maximum p_pi : "<<mod_ppi[1]*P_pi_max<<std::endl;
	std::cout<<"Minimum W : "<<mod_W*W_min<<std::endl;
	std::cout<<"Maximum y : "<<mod_y*y_max<<std::endl;
	std::cout<<"Minimum Mx : "<<mod_Mx*Mx_min<<std::endl;

}

void analyzer::writeCutsToFile(TString outName){
	std::ofstream txtFile;
	txtFile.open(outName);
	for( int i = 0; i < 3; i++ ){
		txtFile<<"Electron Fiducial Region "<<i<<"\t"<<mod_el_fid[i]*e_fid_cuts[i]<<std::endl;
	}
	txtFile<<"PCAL W\t"<<mod_el_PCAL[0]*e_PCAL_W_min<<std::endl;
	txtFile<<"PCAL V\t"<<mod_el_PCAL[1]*e_PCAL_V_min<<std::endl;
	txtFile<<"Electron Min Vertex\t"<<mod_el_vtz[0]*Vz_e_min_inbending<<std::endl;
	txtFile<<"Electron Max Vertex\t"<<mod_el_vtz[1]*Vz_e_max_inbending<<std::endl;
	txtFile<<"Minimum Edep\t"<<mod_el_Edep*e_E_PCAL_min<<std::endl;
	txtFile<<"SF Sigma\t"<<mod_SF_sigma*3.5<<std::endl;	
	txtFile<<"SF Correlation\t"<<mod_el_corr*PCAL_ECIN_SF_min<<std::endl;

	for( int i = 0; i < 3; i++ ){
		txtFile<<"Pip Fiducial Region "<<i<<"\t"<<mod_pi_fid[0][i]*pi_fid_cuts[0][i]<<std::endl;
		txtFile<<"Pim Fiducial Region "<<i<<"\t"<<mod_pi_fid[1][i]*pi_fid_cuts[1][i]<<std::endl;
	}
	txtFile<<"Pip-e Min Vertex\t"<<Vz_pi_mean[mode][0] - 3.5*mod_pi_vtz*Vz_pi_sigma[mode][0]<<std::endl;
	txtFile<<"Pip-e Max Vertex\t"<<Vz_pi_mean[mode][0] + 3.5*mod_pi_vtz*Vz_pi_sigma[mode][0]<<std::endl;

	txtFile<<"Pim-e Min Vertex\t"<<Vz_pi_mean[mode][1] - 3.5*mod_pi_vtz*Vz_pi_sigma[mode][1]<<std::endl;
	txtFile<<"Pim-e Max Vertex\t"<<Vz_pi_mean[mode][1] + 3.5*mod_pi_vtz*Vz_pi_sigma[mode][1]<<std::endl;

	
	
	txtFile<<"Minimum pe\t"<<mod_pe[0]*P_e_min<<std::endl;
	txtFile<<"Maximum pe\t"<<mod_pe[1]*P_e_max<<std::endl;
	txtFile<<"Minimum p_pi\t"<<mod_ppi[0]*P_pi_min<<std::endl;
	txtFile<<"Maximum p_pi\t"<<mod_ppi[1]*P_pi_max<<std::endl;
	txtFile<<"Minimum W\t"<<mod_W*W_min<<std::endl;
	txtFile<<"Maximum y\t"<<mod_y*y_max<<std::endl;
	txtFile<<"Minimum Mx\t"<<mod_Mx*Mx_min<<std::endl;

}

void analyzer::randomizeCuts(){
	std::mt19937 gen(rd());
	
	std::vector<double> cut_vals;
	std::vector<double*> cut_mods;

	for( int i = 0; i < 3; i++ ){
		cut_mods.push_back( &mod_el_fid[i] );
		cut_vals.push_back( e_fid_cuts[i] );
	}

	for( int i = 0; i < 2; i++ ){
		cut_mods.push_back( &mod_el_PCAL[i]);
		cut_mods.push_back( &mod_el_vtz[i]);

		cut_mods.push_back( &mod_pe[i]);
		cut_mods.push_back( &mod_ppi[i]);

		for( int j = 0; j < 3; j++ ){
			cut_mods.push_back( &mod_pi_fid[i][j]);
		}
	}
	
	cut_vals.push_back(e_PCAL_W_min);
	cut_vals.push_back(Vz_e_min_inbending);
	cut_vals.push_back(P_e_min);
	cut_vals.push_back(P_pi_min);
	
	for( int j = 0; j < 3; j++ ){
		cut_vals.push_back( pi_fid_cuts[0][j]);
	}


	cut_vals.push_back(e_PCAL_V_min);
	cut_vals.push_back(Vz_e_max_inbending);
	cut_vals.push_back(P_pi_max);
	cut_vals.push_back(P_e_max);

	for( int j = 0; j < 3; j++ ){
		cut_vals.push_back( pi_fid_cuts[1][j]);
	}
	
	cut_mods.push_back( &mod_el_Edep);
	cut_vals.push_back(e_E_PCAL_min);
	
	cut_mods.push_back( &mod_SF_sigma);
	cut_vals.push_back( 3.5 );
	
	cut_mods.push_back( &mod_el_corr);
	cut_vals.push_back(PCAL_ECIN_SF_min);
	
	cut_mods.push_back( &mod_pi_vtz );
	cut_vals.push_back( 3.5 );

	cut_mods.push_back( &mod_W);
	cut_vals.push_back(W_min);
	
	cut_mods.push_back( &mod_y);
	cut_vals.push_back(y_max);
	
	cut_mods.push_back( &mod_Mx);
	cut_vals.push_back(Mx_min);


	for( int i = 0; i < cut_mods.end() - cut_mods.begin(); i++ ){
		
		double std = cut_vals[i]*0.05;

		std::normal_distribution<double> distribution(cut_vals[i], std);

		*cut_mods[i] = distribution(gen)/(cut_vals[i]);
	}
}

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
		if( e.getEdge(regionIdx) < mod_el_fid[regionIdx]*e_fid_cuts[regionIdx] ){ return false; }
	}	
	return true;
}

bool analyzer::applyElectronPCAL( electron e ){
	return ( e.getW() > mod_el_PCAL[0]*e_PCAL_W_min	
			&&  e.getV() > mod_el_PCAL[1]*e_PCAL_V_min);
}

bool analyzer::applyElectronEDep( electron e ){

	//PCAL MIN EDEP CUT
	return  ( e.getEpcal() > mod_el_Edep*e_E_PCAL_min) ;
}

void analyzer::loadSamplingFractionParams(TString sfFile_name){
	TFile f(sfFile_name);
	for( int sec = 0; sec < 6; sec++ ){
		TF1 * sf_mean = (TF1*)f.Get( Form("fMean_SF_%i", sec) );
		TF1 * sf_sig = (TF1*)f.Get( Form("fSigma_SF_%i", sec) );
		
		for( int par = 0; par < 3; par++ ){
			SF_p_mean[sec][par] = sf_mean->GetParameter(par);
			SF_p_sigma[sec][par] = sf_sig->GetParameter(par);
		}
	}
}

void analyzer::loadSamplingFractionParams(){
	loadSamplingFractionParams( (TString) CUT_PATH + "/SF_fits.root");
}


bool analyzer::applyElectronSF( electron e ){
	//Electron SF cut
	double p = e.get3Momentum().Mag();
	double sf = (e.getEecin() + e.getEpcal() + e.getEecout())/e.get3Momentum().Mag();
	int sec = e.getDC_sector() - 1;
	double mean = SF_p_mean[sec][0] + p*SF_p_mean[sec][1] + p*p*SF_p_mean[sec][2];
	double sigma = SF_p_sigma[sec][0] + p*SF_p_sigma[sec][1] + p*p*SF_p_sigma[sec][2];

	double nSigma = 3.5*mod_SF_sigma;

	bool pass_max = (bool)(sf < ( mean + nSigma*sigma ) );
	bool pass_min = (bool)(sf > ( mean - nSigma*sigma ) );

	return (pass_max && pass_min);
	//return (epid.isElectron(&e)) ;
}
	
bool analyzer::applyElectronCorrelation( electron e ){
	//SF CORRELATION CUT
	if( e.get3Momentum().Mag() < 4.5 )	return true;

	return(  e.getEecin()/e.get3Momentum().Mag() > mod_el_corr*PCAL_ECIN_SF_min - e.getEpcal()/e.get3Momentum().Mag() );
}


bool analyzer::applyElectronVertex( electron e ){
	//ELECTRON VERTEX CUT
	return ((e.getVt().Z() > mod_el_vtz[0]*Vz_e_min_inbending) && (e.getVt().Z() < mod_el_vtz[1]*Vz_e_max_inbending));
}

bool analyzer::applyElectronDetectorCuts( electron e ){
	if (e.getDC_sector() == 0) return false;
	if(!applyElectronFiducials( e ) ) return false;

	// PCAL FIDUCIAL	
	if(!applyElectronPCAL( e ) ) return false;
	//PCAL MIN EDEP CUT
	if(!applyElectronEDep( e ) ) return false;
	//Electron SF cut
	if(!applyElectronSF( e ) ) return false;
	//SF CORRELATION CUT
	if(!applyElectronCorrelation( e ) ) return false;
	//ELECTRON VERTEX CUT
	if(!applyElectronVertex( e ) ) return false;
	
	
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

	

	//double DC_x[3] = {pi.getDC_x1(), pi.getDC_x2(), pi.getDC_x3()};
	//double DC_y[3] = {pi.getDC_y1(), pi.getDC_y2(), pi.getDC_y3()};
	//double DC_z[3] = {pi.getDC_z1(), pi.getDC_z2(), pi.getDC_z3()};

	int piCharge = (int)( pi.getCharge() < 0 );

	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)

		if( pi.getEdge(regionIdx) < mod_pi_fid[piCharge][regionIdx]*pi_fid_cuts[piCharge][regionIdx] ){ 
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
       
	int chargeIdx = (int) ( pi.getCharge() < 0 );
	double min = Vz_pi_mean[mode][chargeIdx] - 3.5*mod_pi_vtz*Vz_pi_sigma[mode][chargeIdx];
	double max = Vz_pi_mean[mode][chargeIdx] + 3.5*mod_pi_vtz*Vz_pi_sigma[mode][chargeIdx];
	
	//DELTA VERTEX CUT
	if( !( (pi.getVt() - e.getVt()).Z() > min && 
				(pi.getVt() - e.getVt()).Z() < max ) ) { 
		return false; 
	}
	
	return true;
}

bool analyzer::applyPionDetectorCuts( pion pi, electron e ){
	
	if (pi.getDC_sector() == 0) return false;
	//Pion DC fiducials
	if(!applyPionDetectorFiducials( pi )) return false;
	//Pion PID Chi2
	if(!applyPionDetectorChi2( pi )) return false;
	//Delta vertex cut
	if(!applyPionDetectorVertex( pi, e )) return false;
       
	return true;
}

//electron kinematical cuts for data and generator
bool analyzer::applyElectronKinematicCuts( electron e ){
		if( sqrt(e.getW2()) < mod_W*W_min ) { return false; }
		if( e.getQ2() < mod_Q2[0]* Q2_min || e.getQ2() > mod_Q2[1]*Q2_max ) { return false; }
                if( e.getXb() < xB_min || e.getXb() > xB_max ) { return false; }
                if( e.getY() > mod_y*y_max ) { return false; }
                
		double p = e.get3Momentum().Mag();
		if( p < mod_pe[0]*P_e_min || p > mod_pe[1]*P_e_max ) { return false; }
                
		double theta = e.get3Momentum().Theta()*rad_to_deg;
		if( theta < theta_min || theta > theta_max ){ return false; }
		
		return true;
}
bool analyzer::applyElectronKinematicCuts( genElectron e ){
		if( sqrt(e.getW2()) < mod_W*W_min ) { return false; }
		if( e.getQ2() < mod_Q2[0]*Q2_min || e.getQ2() > mod_Q2[1]*Q2_max ) { return false; }
                if( e.getXb() < mod_xB[0]*xB_min || e.getXb() > mod_xB[1]*xB_max ) { return false; }
                if( e.getY() > mod_y*y_max ) { return false; }
                
		double p = e.get3Momentum().Mag();
		if( p < mod_pe[0]*P_e_min || p > mod_pe[1]*P_e_max ) { return false; }
                
		double theta = e.get3Momentum().Theta()*rad_to_deg;
		if( theta < theta_min || theta > theta_max ){ return false; }
		
		return true;
}

//Pion kinematical cuts for data and generator
bool analyzer::applyPionKinematicCuts( pion pi ){


	if ( ( pi.getMx() < mod_Mx*Mx_min || pi.getMx() > Mx_max) ) { return false; }
	double p = pi.get3Momentum().Mag();
	if ( p < mod_ppi[0]*P_pi_min || p > mod_ppi[1]*P_pi_max ) { return false; }
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
				match2d[i][j][k]->Print("V");
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
	loadMatchingFunctions3D("matchCut3D.root");
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
	
		if( (sec == 2 || (sec ==3 && p>1)) && phi_temp < 0. ){ phi_temp += 360; }
		cutMin = par_0*(phi_temp - par_1)*(phi_temp - par_1) + par_2;
		
		if( theta > cutMin && theta < 40
			&& theta < max && phi_temp > lower && phi_temp < upper ){ 
			
			return sec; 
		}
	}

	return -1;
}
//Legacy cut function
/*
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
	
	

	//PION CHI2 vs P CUT
	if(! (	( Chi2PID_pion_lowerBound( pi.get3Momentum().Mag(), C ) < pi.getChi2()
         	&& pi.getChi2() < Chi2PID_pion_upperBound( pi.get3Momentum().Mag() , C ) ))) 
		{return false; }
       
	//DELTA VERTEX CUT
	if( !( (pi.getVt() - e.getVt()).Z() > -7 && (pi.getVt() - e.getVt()).Z() < 5 ) ) { return false; }
	//if( ! ( abs( (pi.getVt() - e.getVt()).Z() ) < 20 ) ) { return false; }
	return true;
}
*/
