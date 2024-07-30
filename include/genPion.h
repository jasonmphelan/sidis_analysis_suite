#ifndef __GENPION_H__
#define __GENPION_H__

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "region_particle.h"
#include "rich.h"
#include <vector>

class genPion{// : public TObject {
	public:
		genPion();
		~genPion();


		void Clear();
		void PrintDetectorInfo();
		void PrintKinematicInfo();


		// Particle
		void setSector		(int	iSector		)		{Sector		= iSector	; return;}
		void setPID		(int	iPID		)		{PID		= iPID		; return;}
		void setCharge		(int	iCharge		)		{Charge		= iCharge	; return;}
		
		int	getSector	(void)		{return Sector		;}
		int	getPID		(void)		{return PID		;}
		int	getCharge	(void)		{return Charge		;}

		//Kinematics
		void setMomentum	(clas12::mcpar_ptr mcparts);	
		void setPi_q ( TLorentzVector q, TLorentzVector e4 );
		
		TLorentzVector get4Momentum (void)	{return pi4		;}
		TVector3	get3Momentum (void)	{return pi3		;}	
		TLorentzVector getPi_q (void)		{return pi_q		;}
		TVector3	getVt (void)		{return vt_pi		;}

		void setZ		(double	iZ		)		{Z		= iZ		; return;}
		void setMx		(double iMx		)		{Mx		= iMx		; return;}
		void setXf		(double	iXf		)		{xF		= iXf		; return;}
		void setEta		(double iEta		)		{eta		= iEta		; return;}

		double	getZ		(void)		{return Z		;}
		double	getMx		(void)		{return Mx		;}
		double	getXf		(void)		{return xF		;}
		double	getEta		(void)		{return eta		;}


		// etc
		void setSelection	(bool	iSelection		)	{Selection = iSelection		; return;}
		bool getSelection 	(void)		{return Selection	;}	

		void setKinematicInformation( TLorentzVector q, TLorentzVector pe, clas12::mcpar_ptr rp);

		//ClassDef(clashit,7);

	private:
		int	Sector		;
		int 	PID		;
		int	Charge		;

		TLorentzVector pi4	;
		TVector3 pi3		;
		TLorentzVector pi_q	;
		TVector3 vt_pi		;

		double	Z		;
		double 	Mx		;
		double 	eta		;
		double 	xF		;

		bool	 Selection	;
};



#endif
