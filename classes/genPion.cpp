#include "genPion.h"
#include "constants.h"
#include <iostream>

//ClassImp(clashit);

genPion::genPion(){}	// Empty constructor
genPion::~genPion(){}	// Empty destructor

using namespace constants;

void genPion::Clear(){
	PID		= 0;

	Z		= 0;
	Mx		= 0;
	xF		= 0;
	eta		= 0;

	Selection = false;
	//detectorSelection = false;
	//kinematicSelection = false;
	
	pi4	= TLorentzVector();
	pi3	= TVector3();
	pi_q	= TLorentzVector();
	vt_pi	= TVector3();

}


void genPion::PrintKinematicInfo(){
	double Momentum = pi3.Mag();
	double Theta = pi3.Theta();
	double Phi = pi3.Phi();

	std::cout << "clashit Information Kinematics: Momentum " << Momentum << " , Theta " << Theta << " ,Phi " << Phi << std::endl;
	std::cout << ", Z " << Z << " , Mx " << Mx << std::endl;
}

void genPion::setMomentum ( clas12::mcpar_ptr mcparts){
	pi4.SetXYZM(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), Mpi );

	pi3 = pi4.Vect();
	return;
	
}

void genPion::setPi_q ( TLorentzVector q, TLorentzVector pe ){
    	// move to q-Pe system: q is the z axis, Pe is in x-z plane: Pe=(Pe[x],0,Pe[q])

	TVector3 pi3_temp = pi3;

	pi3_temp.RotateZ( -q.Phi()  );
    	pi3_temp.RotateY( -q.Theta() );
    	pi3_temp.RotateZ( -pe.Phi() );
	
	pi_q.SetVectM( pi3_temp, Mpi ); 
}

void genPion::setKinematicInformation(TLorentzVector q, TLorentzVector pe, clas12::mcpar_ptr rp){
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



