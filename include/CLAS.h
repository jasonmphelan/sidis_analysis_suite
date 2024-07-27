#ifndef CLAS_HH
#define CLAS_HH

class TString;
class TF1;
class TLorentzVector;
class TRandom3;

class CLAS{

public:

	CLAS();
	~CLAS();

	int GetElectronAcceptance(double, double, double);
	void Smear(TLorentzVector *V4, int q); 
	void SmearCD(TLorentzVector *V4, int q); 


private:

	void ReadMomentumParameters(TString, TString);

	void InitMomentumFunctions();
	void InitFiducialFunctions();
	void SetFiducialParameters(double);

	double GetFiducialParameter(int, int, int, double);
	
	double momentumPar[2][6][3][3];	

	TF1* fiducialFunction[2][6];
	TF1* momentumFunction[2][6][3];

	TRandom3 *fRand;
};

#endif

