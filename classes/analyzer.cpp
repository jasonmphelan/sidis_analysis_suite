// Erez O. C., Oct-6, 2021
#include "analyzer.h"
#define CUT_PATH _DATA


analyzer::analyzer(int _fdebug_, int _torusBending_){
    SetVerbosity    (_fdebug_);
    SetTorusBending (_torusBending_);
}

analyzer::~analyzer(){}

bool 	applyElectronDetectorCuts( electron e );
//bool 	applyElectronKinematics( electron e );

bool 	applyPionKinematics( pion pi );
//bool 	applyPionDetectorCuts( pion pi );

//void	loadAcceptanceMatching (int torusBending=1); //  -1 for In-bending, +1 for Out-bending
//void	loadCorrections (TFile corrFileName); //  -1 for In-bending, +1 for Out-bending
//void	getCorrections ( electron e, pion pi);
//void	applyAcceptanceMatching( pion pi );	
void	printCutValues ();
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



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
double analyzer::ComputeLightConeFraction( TLorentzVector p ){
    // compute light-cone momentum fraction
    double m = p.Mag();
    double alpha = (p.E() - p.Z())/m;
    return alpha;
}


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



