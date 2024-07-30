#include "genElectron.h"
#include "constants.h"
#include <iostream>

//ClassImp(clashit);

genElectron::genElectron(){}	// Empty constructor
genElectron::~genElectron(){}	// Empty destructor

using namespace constants;

void genElectron::Clear(){
	Sector		= -1;
	PID		= 0;
	Charge		= 0;

	Q2		= 0;
	Omega		= 0;
	Xb		= 0;
	W2		= 0;
	y		= 0;

	Selection = false;
	//detectorSelection = false;
	//kinematicSelection = false;
	
	e4	= TLorentzVector();
	e3	= TVector3();
	q	= TLorentzVector();
	vt	= TVector3();

}

void genElectron::PrintKinematicInfo(){
	double Momentum = e3.Mag();
	double Theta = e3.Theta();
	double Phi = e3.Phi();

	std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi;
	std::cout << ", Q2 " << Q2 << " , Omega/nu " << Omega << " , Xb " << Xb << " , W2 " << W2 << std::endl;
	//std::cout << "clashit Information q-vector: Magitude(q) " << Q << " , ThetaQ " << ThetaQ << " ,PhiQ " << PhiQ << std::endl;
	
}


void genElectron::setMomentum ( clas12::mcpar_ptr mcparts){
	e4.SetXYZM(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), Me );

	e3 = e4.Vect();
	return;
	
}

void genElectron::setKinematicInformation( double Ebeam, clas12::mcpar_ptr rp){
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );
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



