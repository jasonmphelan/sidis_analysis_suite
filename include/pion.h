#ifndef __ELECTRON_H__
#define __ELECTRON_H__

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "region_particle.h"
#include <vector>

class pion{// : public TObject {
	public:
		pion();
		~pion();


		void Clear();
		void Print();


		// Particle
		void setSector		(int	iSector		)		{Sector		= iSector	; return;}
		void setPID		(int	iPID		)		{PID		= iPID		; return;}
		void setCharge		(int	iCharge		)		{Charge		= iCharge	; return;}
		void setStatus		(int	iStatus		)		{Status		= iStatus	; return;}
		
		int	getSector	(void)		{return Sector		;}
		int	getPID		(void)		{return PID		;}
		int	getCharge	(void)		{return Charge		;}
		int	getStatus	(void)		{return Status		;}

		void setTime		(double	iTime		)		{Time		= iTime		; return;}
		void setBeta		(double	iBeta		)		{Beta		= iBeta		; return;}
		void setChi2		(double	iChi2		)		{Chi2		= iChi2		; return;}
		
		double	getTime		(void)		{return Time		;}
		double	getBeta		(void)		{return Beta		;}
		double	getChi2		(void)		{return Chi2		;}
		
		// Calorimeter
		void setEtot		(double	iEtot		)		{Etot		= iEtot		; return;}
		void setEpcal		(double	iEpcal		)		{Epcal		= iEpcal	; return;}
		void setEecin		(double	iEecin		)		{Eecin		= iEecin	; return;}
		void setEecout		(double	iEecout		)		{Eecout		= iEecout	; return;}
		void setEoP		(double	iEoP		)		{EoP		= iEoP		; return;}
		void setTimeScint	(double	iTimeScint	)		{TimeScint	= iTimeScint	; return;}
		void setPathScint	(double	iPathScint	)		{PathScint	= iPathScint	; return;}
		void setU		(double	iU		)		{U		= iU		; return;}
		void setV		(double	iV		)		{V		= iV		; return;}
		void setW		(double	iW		)		{W		= iW		; return;}
		void setPCal_X		(double	iX		)		{PCal_X		= iX		; return;}
		void setPCal_Y		(double	iY		)		{PCal_Y		= iY		; return;}
		void setPCal_Z		(double	iZ		)		{PCal_Z		= iZ		; return;}
		
		double	getEtot		(void)		{return Etot		;}
		double	getEpcal	(void)		{return Epcal		;}
		double	getEecin	(void)		{return Eecin		;}
		double	getEecout	(void)		{return Eecout		;}
		double	getEoP		(void)		{return EoP		;}
		double	getTimeScint	(void)		{return TimeScint	;}
		double	getPathScint	(void)		{return PathScint	;}
		double	getU		(void)		{return U		;}
		double	getV		(void)		{return V		;}
		double	getW		(void)		{return W		;}
		double	getPCal_X	(void)		{return PCal_X		;}
		double	getPCal_Y	(void)		{return PCal_Y		;}
		double	getPCal_Z	(void)		{return PCal_Z		;}
		
		//Kinematics
		void setMomentum	(clas12::region_part_ptr rp);	
		void setVt		(clas12::region_part_ptr rp);
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


		// Drift Chamber
		void setDC_chi2   	(double iDC_chi2   	)		{DC_chi2	= iDC_chi2	; return;}
		void setDC_NDF		(int    iDC_NDF		)		{DC_NDF 	= iDC_NDF 	; return;}
		void setDC_sector 	(int    iDC_sector  	) 		{DC_sector      = iDC_sector    ; return;}

		void setDC_x1		(double iDC_x1		)		{DC_x1   	= iDC_x1 	; return;}
		void setDC_y1		(double iDC_y1		)		{DC_y1   	= iDC_y1   	; return;}
		void setDC_z1		(double iDC_z1		)		{DC_z1   	= iDC_z1   	; return;}

		void setDC_x2		(double iDC_x2		)		{DC_x2   	= iDC_x2 	; return;}
		void setDC_y2		(double iDC_y2		)		{DC_y2   	= iDC_y2   	; return;}
		void setDC_z2		(double iDC_z2		)		{DC_z2   	= iDC_z2   	; return;}

		void setDC_x3		(double iDC_x3		)		{DC_x3   	= iDC_x3 	; return;}
		void setDC_y3		(double iDC_y3		)		{DC_y3   	= iDC_y3   	; return;}
		void setDC_z3		(double iDC_z3		)		{DC_z3   	= iDC_z3   	; return;}
		
		double  getDC_chi2      (void)          {return DC_chi2         ;}
		int     getDC_NDF       (void)          {return DC_NDF          ;}
		int     getDC_sector    (void)          {return DC_sector       ;}

		double  getDC_x1        (void)          {return DC_x1           ;}
		double  getDC_y1        (void)          {return DC_y1           ;}
		double  getDC_z1        (void)          {return DC_z1           ;}

		double  getDC_x2        (void)          {return DC_x2           ;}
		double  getDC_y2        (void)          {return DC_y2           ;}
		double  getDC_z2        (void)          {return DC_z2           ;}

		double  getDC_x3        (void)          {return DC_x3           ;}
		double  getDC_y3        (void)          {return DC_y3           ;}
		double  getDC_z3        (void)          {return DC_z3           ;}
	
		// Cherenkov banks
		void setNphe		(double iNphe		)		{Nphe		= iNphe		; return;}
		void setKov_x		(double iKov_x		)		{Kov_x		= iKov_x	; return;}
		void setKov_y		(double iKov_y		)		{Kov_y		= iKov_y	; return;}
		void setKov_z		(double iKov_z		)		{Kov_z		= iKov_z	; return;}
		void setKov_chi2	(double iKov_chi2	)		{Kov_chi2	= iKov_chi2	; return;}
		void setKov_time	(double iKov_time	)		{Kov_time	= iKov_time	; return;}
		void setKov_path	(double iKov_path	)		{Kov_path	= iKov_path	; return;}
		void setKov_det		(int    iKov_det	)		{Kov_det	= iKov_det	; return;}
		void setKov_sector	(int    iKov_sec	)		{Kov_sec	= iKov_sec	; return;}
		void setKov_status	(int    iKov_status	)		{Kov_status	= iKov_status	; return;}
		
		double getNphe		(void)		{return Nphe		;}
		double getKov_x		(void)		{return Kov_x		;}
		double getKov_y		(void)		{return Kov_y		;}
		double getKov_z		(void)		{return Kov_z		;}
		double getKov_chi2	(void)		{return Kov_chi2	;}
		double getKov_time	(void)		{return Kov_time	;}
		double getKov_path	(void)		{return Kov_path	;}
		int    getKov_det	(void)		{return Kov_det		;}
		int    getKov_sector	(void)		{return Kov_sec		;}
		int    getKov_status	(void)		{return Kov_status	;}
	
		// Scintillator banks
		void setScint_status	(int    iScint_status		)	{Scint_status		.push_back( iScint_status	)	; return; }
		void setScint_sector	(int    iScint_sector		)	{Scint_sector		.push_back( iScint_sector	)	; return; }
		void setScint_layer	(int    iScint_layer		)	{Scint_layer		.push_back( iScint_layer	)	; return; }
		void setScint_component	(int    iScint_component	)	{Scint_component	.push_back( iScint_component	)	; return; }
		void setScint_Edep	(double iScint_Edep		)	{Scint_Edep		.push_back( iScint_Edep		)	; return; }
		void setScint_time	(double iScint_time		)	{Scint_time		.push_back( iScint_time		)	; return; }
		void setScint_path	(double iScint_path		)	{Scint_path		.push_back( iScint_path		)	; return; }
		void setScint_chi2	(double iScint_chi2		)	{Scint_chi2		.push_back( iScint_chi2		)	; return; }
		void setScint_x		(double iScint_x		)	{Scint_x		.push_back( iScint_x		)	; return; }
		void setScint_y		(double iScint_y		)	{Scint_y		.push_back( iScint_y		)	; return; }
		void setScint_z		(double iScint_z		)	{Scint_z		.push_back( iScint_z		)	; return; }
		
		std::vector<int>    getScint_status		(void)		{return Scint_status	;}
		std::vector<int>    getScint_sector		(void)		{return Scint_sector	;}
		std::vector<int>    getScint_layer		(void)		{return Scint_layer	;}
		std::vector<int>    getScint_component		(void)		{return Scint_component	;}
		std::vector<double> getScint_Edep		(void)		{return Scint_Edep	;}
		std::vector<double> getScint_time		(void)		{return Scint_time	;}
		std::vector<double> getScint_path		(void)		{return Scint_path	;}
		std::vector<double> getScint_chi2		(void)		{return Scint_chi2	;}
		std::vector<double> getScint_x			(void)		{return Scint_x		;}
		std::vector<double> getScint_y			(void)		{return Scint_y		;}
		std::vector<double> getScint_z			(void)		{return Scint_z		;}


		//Rich
		void setBeta_rich 	(double iBeta			)	{Beta_rich = iBeta		; return;}

		double getBeta_rich	(void)		{return Beta_rich	;}
		
		// etc
		void setSelection	(bool	iSelection		)	{Selection = iSelection		; return;}
		bool getSelection 	(void)		{return Selection	;}	

		void setDetectorInformation	(clas12::region_part_ptr);
		void setKinematicInformation( TLorentzVector q, TLorentzVector pe, clas12::region_part_ptr rp);
		void setPion( TLorentzVector q, TLorentzVector pe, clas12::region_part_ptr rp );

		//ClassDef(clashit,7);

	private:
		int	Sector		;
		int 	PID		;
		int	Charge		;
		int	Status		;

		double	Time		;
		double	Beta		;
		double	Chi2		;
		double	Etot		;
		double	Epcal		;
		double	Eecin		;
		double	Eecout		;
		double	EoP		;
		double	TimeScint	;
		double	PathScint	;
		double 	U		;
		double	V		;
		double  W		;
		double 	PCal_X		;
		double	PCal_Y		;
		double  PCal_Z		;

		TLorentzVector pi4	;
		TVector3 pi3		;
		TLorentzVector pi_q	;
		TVector3 vt_pi		;

		double	Z		;
		double 	Mx		;
		double 	eta		;
		double 	xF		;

		double DC_chi2          ;
		int    DC_NDF           ;
		int    DC_sector        ;

		double DC_x1            ;
		double DC_y1            ;
		double DC_z1            ;

		double DC_x2            ;
		double DC_y2            ;
		double DC_z2            ;

		double DC_x3            ;
		double DC_y3            ;
		double DC_z3            ;

		double Nphe		;
                double Kov_x		;
                double Kov_y		;
                double Kov_z		;
                double Kov_chi2		;
                double Kov_time		;
                double Kov_path		;
                int    Kov_det		;
                int    Kov_sec		;
                int    Kov_status	;

		std::vector<int>	Scint_status	;
		std::vector<int>	Scint_sector	;
		std::vector<int>	Scint_layer	;
		std::vector<int>	Scint_component	;
		std::vector<double>	Scint_Edep	;
		std::vector<double>	Scint_time	;
		std::vector<double>	Scint_path	;
		std::vector<double>	Scint_chi2	;
		std::vector<double>	Scint_x		;
		std::vector<double>	Scint_y		;
		std::vector<double>	Scint_z		;

		bool	 Selection	;
		double	Beta_rich	;
};



#endif
