#ifndef __GENELECTRON_H__
#define __GENELECTRON_H__

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "region_particle.h"
#include <vector>
#include "mcparticle.h"

class genElectron{// : public TObject {
	public:
		genElectron();
		~genElectron();


		void Clear();
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
		void setQ		(TLorentzVector iQ){ q = iQ; };
		
		TLorentzVector get4Momentum (void)	{return e4		;}
		TVector3	get3Momentum (void)	{return e3		;}	
		TLorentzVector getQ (void)		{return q		;}

		void setQ2		(double	iQ2		)		{Q2		= iQ2		; return;}
		void setOmega		(double	iOmega		)		{Omega		= iOmega	; return;}
		void setXb		(double	iXb		)		{Xb		= iXb		; return;}
		void setW2		(double	iW2		)		{W2		= iW2		; return;}
		void setY		(double iY		)		{y		= iY		; return;}

		double	getQ2		(void)		{return Q2		;}
		double	getOmega	(void)		{return Omega		;}
		double	getXb		(void)		{return Xb		;}
		double	getW2		(void)		{return W2		;}
		double	getY		(void)		{return y		;}


		
		// etc
		void setSelection	(bool	iSelection		)	{Selection = iSelection		; return;}
		bool getSelection 	(void)		{return Selection	;}	

		void setKinematicInformation( double Ebeam, clas12::mcpar_ptr rp);

		//ClassDef(clashit,7);

	private:
		int	Sector		;
		int 	PID		;
		int	Charge		;

		TLorentzVector e4	;
		TVector3 e3		;
		TLorentzVector q	;
		TVector3 vt		;

		double 	Q2		;
		double 	Omega		;
		double	Xb		;
		double	W2		;
		double 	y		;


		bool	 Selection	;
};



#endif
