#include "electron.h"
#include "constants.h"
#include <iostream>

//ClassImp(clashit);

electron::electron(){}	// Empty constructor
electron::~electron(){}	// Empty destructor

using namespace constants;

void electron::Clear(){
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

	Q2		= 0;
	Omega		= 0;
	Xb		= 0;
	W2		= 0;
	y		= 0;

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
/*
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
*/
	Selection = false;
	//detectorSelection = false;
	//kinematicSelection = false;
	
	e4	= TLorentzVector();
	e3	= TVector3();
	q	= TLorentzVector();
	vt	= TVector3();

	Beta_rich	= 0;
}

void electron::PrintDetectorInfo(){

	std::cout << "clashit Information REC: PID " << PID << " , Charge " << Charge << " , Status " << Status;
	std::cout << ", Sector(Calo) " << Sector << " , Chi2 " << Chi2 << ", Time " << Time << " , Beta " << Beta << std::endl;
	std::cout << "clashit Information CALO: Etot " << Etot << " , Epcal " << Epcal << ", Eecin " << Eecin;
	std::cout << ", Eecout " << Eecout << ", EoverP " << EoP << " , U(PCal) " << U << " , V(ECal) " << V << " , W(ECal) " << W << std::endl;
	//std::cout << "clashit Information Vertex: Vtx " << Vtx << " , Vty " << Vty << ", Vtz " << Vtz;
	//std::cout << ", TimeScint(-starttime) " << TimeScint << " , Pathlength Scint " << PathScint << std::endl;
	std::cout << "clashit Information DC: DC_chi2 " << DC_chi2 << " , DC_NDF " << DC_NDF << ", DC_sector " << DC_sector;
	std::cout << ", DC_x1 " << DC_x1 << " , DC_y1 " << DC_y1 << " , DC_z1 " << DC_z1 << std::endl;
	std::cout << "clashit Information DC: DC_x2 " << DC_x2 << " , DC_y2 " << DC_y2<< " , DC_z2 " << DC_z2;
	std::cout << ", DC_x3 " << DC_x3 << " , DC_y3 " << DC_y3 << " , DC_z3 " << DC_z3 << std::endl;
}

void electron::PrintKinematicInfo(){
	double Momentum = e3.Mag();
	double Theta = e3.Theta();
	double Phi = e3.Phi();

	std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi;
	std::cout << ", Q2 " << Q2 << " , Omega/nu " << Omega << " , Xb " << Xb << " , W2 " << W2 << std::endl;
	//std::cout << "clashit Information q-vector: Magitude(q) " << Q << " , ThetaQ " << ThetaQ << " ,PhiQ " << PhiQ << std::endl;
	
}

void electron::setVt(clas12::region_part_ptr rp){
	TVector3 V(rp->par()->getVx(),
			rp->par()->getVy(),
			rp->par()->getVz());
	vt = V;
	return;
}

void electron::setMomentum ( clas12::region_part_ptr rp){

	e4.SetXYZM(rp->par()->getPx(),
			rp->par()->getPy(),
			rp->par()->getPz(),
			Me);
	e3 = e4.Vect();
	return;
	
}

void electron::setKinematicInformation( double Ebeam, clas12::region_part_ptr rp){
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );
	setVt(rp);
	setMomentum(rp);

    	setQ( beam - e4 );
    	setQ2( (double) -q.Mag2() );
	setOmega( q.E() );
    	setXb( Q2/(2. * Mp * q.E()) );
	setY( Omega / Ebeam );
    	setW2( (double)((p_rest + q).Mag2()) );
    	//setW_d( sqrt((d_rest + q).Mag2()) );
	return;
}


void electron::setDetectorInformation(clas12::region_part_ptr rp){
	//setCharge( rp->par()->getCharge() );

	//check if momentum set:
	if( e4.P() == 0 ){ setMomentum(rp); }

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
	setEoP((getEpcal() + getEecin() + getEecout())/e4.P());

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

	/*
	// RICH info
	auto RICH_info = rp->rich();
	double temp_rich_angle = RICH_info->getBest_ch();

	if( temp_rich_angle > 0. ){
		e.setRich_beta( 1./(rich_n*cos(temp_rich_angle)) );
	}
	else{ e.setRich_beta(0.); }
	*/
	return;
}

void electron::setElectron( double Ebeam, clas12::region_part_ptr rp ){
	setKinematicInformation(Ebeam, rp);
	setDetectorInformation(rp);
	return;
}

