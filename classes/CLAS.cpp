#include "CLAS.h"

#include "TF1.h"
#include "TString.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include <fstream>
#include <unistd.h>

#define DAT_PATH _DATA

using namespace std;

CLAS::CLAS() {

	ReadMomentumParameters((TString) DAT_PATH+"/acceptance_map/upper_momentum_fit.dat", (TString) DAT_PATH+"/acceptance_map/lower_momentum_fit.dat");
	//ReadMomentumParameters("./dat/upper_momentum_fit.dat", "./dat/lower_momentum_fit.dat");
	InitMomentumFunctions();
	InitFiducialFunctions();

	fRand = new TRandom3(0);

}

CLAS::~CLAS() {
}

void CLAS::ReadMomentumParameters(TString upperParFile, TString lowerParFile) {

	ifstream upperFile(upperParFile);
	ifstream lowerFile(lowerParFile);
	
	int sector, fidPar;
	double momPar0, momPar1, momPar2;

	while(upperFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[0][sector-1][fidPar][0] = momPar0;
		momentumPar[0][sector-1][fidPar][1] = momPar1;
		momentumPar[0][sector-1][fidPar][2] = momPar2;

	}

	while(lowerFile >> sector >> fidPar >> momPar0 >> momPar1 >> momPar2) {

		momentumPar[1][sector-1][fidPar][0] = momPar0;
		momentumPar[1][sector-1][fidPar][1] = momPar1;
		momentumPar[1][sector-1][fidPar][2] = momPar2;

	}


}

void CLAS::InitMomentumFunctions() {

	TString momentum_functions[3] = {"[0] + exp([1]*x+[2])", "[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int fidPar = 0; fidPar < 3; fidPar++) {
			for(int i = 0; i < 2; i++) {
				momentumFunction[i][sector][fidPar] = new TF1(Form("%sMomFunc_%i_%i", boundary[i].Data(), sector, fidPar), momentum_functions[fidPar], 0., 12.);
				for(int momPar = 0; momPar < 3; momPar++) {
					momentumFunction[i][sector][fidPar]->SetParameter(momPar, momentumPar[i][sector][fidPar][momPar]);
				}	
			}
		}
	}
}


void CLAS::InitFiducialFunctions() {

	TString fiducial_functions[2] = {"[0] - exp([1]*x+[2])", "[0] + exp([1]*x+[2])"};
	TString boundary[2] = {"upper","lower"};

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			fiducialFunction[i][sector] = new TF1(Form("%sFidFunc_%i", boundary[i].Data(), sector), fiducial_functions[i], 0., 60.);
		}
	}
}

void CLAS::SetFiducialParameters(double pe) {

	for(int sector = 0; sector < 6; sector++) {
		for(int i = 0; i < 2; i++) {
			for(int fidPar = 0; fidPar < 3; fidPar++) {
				double par = GetFiducialParameter(i, sector, fidPar, pe);
				fiducialFunction[i][sector]->SetParameter(fidPar, par);
			}
		}
	}
}


double CLAS::GetFiducialParameter(int boundary, int sector, int par, double pe) {

	return momentumFunction[boundary][sector][par]->Eval(pe);

}


int CLAS::GetElectronAcceptance(double theta, double phi, double pe) {

	SetFiducialParameters(pe);

	for(int sector = 0; sector < 6; sector++) {

		double maxPhi = fiducialFunction[0][sector]->Eval(theta);
		double minPhi = fiducialFunction[1][sector]->Eval(theta);

		if (phi > minPhi && phi < maxPhi) {
			return (sector + 1);
		}
	}

	// if it's not already in a sector, 
	// check for values of (phi - 360) 
	 for(int sector = 0; sector < 6; sector++) {

                double maxPhi = fiducialFunction[0][sector]->Eval(theta);
                double minPhi = fiducialFunction[1][sector]->Eval(theta);

                if ((phi - 360.) > minPhi && (phi - 360.) < maxPhi) {
                        return (sector + 1);
                }
        }


	return -1;
}


void CLAS::Smear(TLorentzVector *V4, int q){
        double inM = V4->M();
        double sP  = V4->P();
        double sTh = V4->Theta();
        double sPh = V4->Phi();

        double sThD = TMath::RadToDeg()*sTh;
        double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
        double momS2 = 0.02*sThD;
        double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
        momR *= 2.0;

        double theS1 = 0.004*sThD + 0.1;
        double theS2 = 0;
        double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
        theR *= 2.5;

        double phiS1 = 0.85-0.015*sThD;
        double phiS2 = 0.17-0.003*sThD;
        double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
        phiR *= 3.5;

        sPh += TMath::DegToRad() * phiR * fRand->Gaus();
        sTh += TMath::DegToRad() * theR * fRand->Gaus();
        sP  += momR  * fRand->Gaus() *  V4->P() ;

        V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
        V4->SetRho( sP );
        V4->SetTheta( sTh );
        V4->SetPhi( sPh );
}

void CLAS::SmearCD(TLorentzVector *V4, int q){
        double inM = V4->M();
        double sP  = V4->P();
        double sTh = V4->Theta();
        double sPh = V4->Phi();

        double sThD = TMath::RadToDeg()*sTh;
        double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
        double momS2 = 0.02*sThD;
        double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
        momR *= 5;

        double theS1 = 0.004*sThD + 0.1;
        double theS2 = 0;
        double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
        theR *= 1.0;

        double phiS1 = 0.85-0.015*sThD;
        double phiS2 = 0.17-0.003*sThD;
        double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
        phiR *= 3.5;

        sPh += TMath::DegToRad() * phiR * fRand->Gaus();
        sTh += TMath::DegToRad() * theR * fRand->Gaus();
        sP  += momR  * fRand->Gaus() *  V4->P() ;

        V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
        V4->SetRho( sP );
        V4->SetTheta( sTh );
        V4->SetPhi( sPh );
}
