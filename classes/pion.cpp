#include "pion.h"
#include "constants.h"
#include <iostream>

//ClassImp(clashit);

pion::pion(){}	// Empty constructor
pion::~pion(){}	// Empty destructor

using namespace constants;

void pion::Clear(){
	Sector		= -1;
	PID		= 0;
	Charge		= 0;
	Status		= 0;

	Time		= 0;
	Beta		= 0;
	Chi2		= 0;
	Etot		= 0;
	Epcal		= 0;
	Eecin		= 0;
	Eecout		= 0;
	EoP		= 0;
	TimeScint	= 0;
	PathScint	= 0;
	U		= 0;
	V		= 0;
	W		= 0;
	PCal_X		= 0;
	PCal_Y		= 0;
	PCal_Z		= 0;

	Z		= 0;
	Mx		= 0;
	xF		= 0;
	eta		= 0;

	DC_chi2         = -999;
	DC_NDF          = -999;
	DC_sector       = -999;

	DC_x1           = -999;
	DC_y1           = -999;
	DC_z1           = -999;

	DC_x2           = -999;
	DC_y2           = -999;
	DC_z2           = -999;

	DC_x3           = -999;
	DC_y3           = -999;
	DC_z3           = -999;

	Nphe		= -999;
	Kov_x		= -999;
	Kov_y		= -999;
	Kov_z		= -999;
	Kov_chi2	= -999;
	Kov_time	= -999;
	Kov_path	= -999;
	Kov_det		= -999;
	Kov_sec		= -999;
	Kov_status	= -999;

	Scint_status	.clear();
	Scint_sector	.clear();
	Scint_layer	.clear();
	Scint_component	.clear();
	Scint_Edep	.clear();
	Scint_time	.clear();
	Scint_path	.clear();
	Scint_chi2	.clear();
	Scint_x		.clear();
	Scint_y		.clear();
	Scint_z		.clear();

	Selection = false;
	//detectorSelection = false;
	//kinematicSelection = false;
	
	pi4	= TLorentzVector();
	pi3	= TVector3();
	pi_q	= TLorentzVector();
	vt_pi	= TVector3();

	Beta_rich	= 0;
}

void pion::PrintDetectorInfo(){

	std::cout << "clashit Information REC: PID " << PID << " , Charge " << Charge << " , Status " << Status;
	std::cout << ", Sector(Calo) " << Sector << " , Chi2 " << Chi2 << ", Time " << Time << " , Beta " << Beta << std::endl;
	std::cout << "clashit Information CALO: Etot " << Etot << " , Epcal " << Epcal << ", Eecin " << Eecin;
	std::cout << ", Eecout " << Eecout << ", EoverP " << EoP << " , U(PCal) " << U << " , V(ECal) " << V << " , W(ECal) " << W << std::endl;
	//std::cout << "clashit Information Vertex: Vtx " << Vtx << " , Vty " << Vty << ", Vtz " << Vtz;
	std::cout << ", TimeScint(-starttime) " << TimeScint << " , Pathlength Scint " << PathScint << std::endl;
	std::cout << "clashit Information DC: DC_chi2 " << DC_chi2 << " , DC_NDF " << DC_NDF << ", DC_sector " << DC_sector;
	std::cout << ", DC_x1 " << DC_x1 << " , DC_y1 " << DC_y1 << " , DC_z1 " << DC_z1 << std::endl;
	std::cout << "clashit Information DC: DC_x2 " << DC_x2 << " , DC_y2 " << DC_y2<< " , DC_z2 " << DC_z2;
	std::cout << ", DC_x3 " << DC_x3 << " , DC_y3 " << DC_y3 << " , DC_z3 " << DC_z3 << std::endl;
	std::cout << "RICH information: Beta "<<Beta_rich<<std::endl;
	//std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi;
	//std::cout << ", Q2 " << Q2 << " , Omega/nu " << Omega << " , Xb " << Xb << " , W2 " << W2 << std::endl;
	//std::cout << "clashit Information q-vector: Magitude(q) " << Q << " , ThetaQ " << ThetaQ << " ,PhiQ " << PhiQ << std::endl;
	
}

void pion::PrintKinematicInfo(){
	double Momentum = pi3.Mag();
	double Theta = pi3.Theta();
	double Phi = pi3.Phi();

	std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi << std::endl;
	std::cout << ", Z " << Z << " , Mx " << Mx << std::endl;
}

void pion::setVt(clas12::region_part_ptr rp){
	TVector3 V(rp->par()->getVx(),
			rp->par()->getVy(),
			rp->par()->getVz());
	vt_pi = V;
	return;
}

void pion::setMomentum ( clas12::region_part_ptr rp){

	pi4.SetXYZM(rp->par()->getPx(),
			rp->par()->getPy(),
			rp->par()->getPz(),
			Mpi);
	pi3 = pi4.Vect();
	return;
	
}

void pion::setPi_q ( TLorentzVector q, TLorentzVector pe ){
    	// move to q-Pe system: q is the z axis, Pe is in x-z plane: Pe=(Pe[x],0,Pe[q])

	TVector3 pi3_temp = pi3;

	pi3_temp.RotateZ( -q.Phi()  );
    	pi3_temp.RotateY( -q.Theta() );
    	pi3_temp.RotateZ( -pe.Phi() );
	
	pi_q.SetVectM( pi3_temp, Mpi ); 
}

void pion::setKinematicInformation(TLorentzVector q, TLorentzVector pe, clas12::region_part_ptr rp){
	setVt(rp);
	
	if( pi4.P() == 0 ){ setMomentum(rp); }

	setPi_q( q, pe ); 					

	setZ(		pi4.E()/q.E()	);
	//Z_LC.push_back(		(pi_q_dummy.E() + pi_q_dummy.Pz()) / (q.E() + q.P())	);

	setMx(		( q + p_rest - pi4 ).Mag()	);
	//M_x_d.push_back(	( q + d_rest - p_dummy ).Mag()	);

	//setXf(		2. * (pi3.Dot(q.Vect())) / ( q.P() * sqrt(e.getW2) )	);

	setEta(		0.5 * log((pi_q.E()+pi_q.Pz()) /
						(pi_q.E()-pi_q.Pz()))	);
}


void pion::setDetectorInformation(clas12::region_part_ptr rp){
	setCharge( rp->par()->getCharge() );
	setPID ( (int) (211*Charge) );

	//check if momentum set:
	if( pi4.P() == 0 ){ setMomentum(rp); }

	setChi2(rp->par()->getChi2Pid());

	// detector information on electron
	auto e_PCAL_info= rp->cal(clas12::PCAL);
	setEpcal(e_PCAL_info->getEnergy());
	setSector(e_PCAL_info->getSector());
	setV( e_PCAL_info->getLv());
	setW(e_PCAL_info->getLw());
	setEecin(rp->cal(clas12::ECIN)->getEnergy());
	setEecout(rp->cal(clas12::ECOUT)->getEnergy());

	// hit position in PCAL
	setPCal_X(e_PCAL_info->getX());
	setPCal_Y( e_PCAL_info->getY());
	setPCal_Z(e_PCAL_info->getZ());

	// Sampling Fraction		
	setEoP((getEpcal() + getEecin() + getEecout())/pi4.P());

	// Drift Chamber tracking system
	auto e_DC_info  = rp->trk(clas12::DC);
	setDC_sector(e_DC_info->getSector()); // tracking sector
	setDC_chi2(e_DC_info->getChi2());  // tracking chi^2/NDF
	setDC_NDF(e_DC_info->getNDF());  // tracking chi^2/NDF

	setDC_x1(rp->traj(clas12::DC,DC_layers[0])->getX());
	setDC_y1(rp->traj(clas12::DC,DC_layers[0])->getY());	
	setDC_z1(rp->traj(clas12::DC,DC_layers[0])->getZ());

	setDC_x2(rp->traj(clas12::DC,DC_layers[1])->getX());
	setDC_y2(rp->traj(clas12::DC,DC_layers[1])->getY());	
	setDC_z2(rp->traj(clas12::DC,DC_layers[1])->getZ());

	setDC_x3(rp->traj(clas12::DC,DC_layers[2])->getX());
	setDC_y3(rp->traj(clas12::DC,DC_layers[2])->getY());	
	setDC_z3(rp->traj(clas12::DC,DC_layers[2])->getZ());
	setBeta(rp->par()->getBeta());

	
	// RICH info
	auto RICH_info = rp->rich();
	double temp_rich_angle = RICH_info->getBest_ch();

	if( temp_rich_angle > 0. ){
		setBeta_rich( 1./(n_rich*cos(temp_rich_angle)) );
	}
	else{ setBeta_rich(0.); }
	
	return;
}

void pion::setPion( TLorentzVector q, TLorentzVector pe,  clas12::region_part_ptr rp ){
	setKinematicInformation(q, pe , rp);
	setDetectorInformation(rp);
	return;
}

