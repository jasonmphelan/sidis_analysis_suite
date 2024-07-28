#ifndef READER_HH
#define READER_HH

#include "HipoChain.h"
#include "TString.h"

class reader{

public:

	reader();
	~reader();
	
	void setNumFiles( int nFiles_i ){ nFiles = nFiles_i; }
	int getRunNum( int i );
	void setRunType( int i ){ runType = i; }
	void setEnergy( double i ){ EBeam = i; }
	void readRunFiles( clas12root::HipoChain &fileList);
	
	void getRunFiles( clas12root::HipoChain &files);
	void getRunSkimsByName( TChain * chain, TString name );
	void getRunSkimsAllEnergy( TChain * chain, TString name );
	void readRunFilesAllE(clas12root::HipoChain &fileList);

private:

	void getSkimsByName( TChain * chain, TString name );
	
	void setDataPaths();

	TString dataPath;
	TString runList;

	std::vector<int> runNums;
	void getRunList();

	int nFiles;
	int runType = 0;
	double EBeam = 10.2;
};

#endif

