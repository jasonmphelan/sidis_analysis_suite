#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
#include "reader.h"
#include "analyzer.h"
#include "reader.h"
#include "TArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace clas12;
using namespace constants;

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
Cuts applied in order
1) PCAL FID
2) PCAL MIN
3) SF
4) SF Correlation
5) Vz_e
6) Electron DC Fiducials
7) Vz_pi
8) Chi2
9) DC Fiducials
*/

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [outFile name] [Beam energy (0 for all energies)]\n";
		return -1;
	}
	
       	double Ebeam = atof(argv[2]); // [GeV]
	TString outFile_name = argv[1]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 && Ebeam != 0 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	anal.setAnalyzerLevel(0);
	anal.loadCutValues(-1, Ebeam);
	
	reader runReader;
	runReader.setNumFiles( 5 );
	runReader.setRunType( 0 );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
       	runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;

	//Declare hists
	TString spec[2] = {"pip", "pim"};

	TH2F * hSF[6];
	TF1 * fMean[6];
	TF1 * fSigma[6];
	
	TGraph *  gThetaPMax[6][2];
	TGraph *  gThetaPMin[6][2];

	for( int j = 0; j < 6; j++ ){
		hSF[j] = new TH2F( Form("hSF_sec_%i", j), ";p_{e} [GeV];#Delta(E_{PCAL} + E_{ECIN} + E_{ECOUT})/p_{e}", 150, 0, 8, 150, .1, .35 );
		fMean[j] = new TF1( Form("fMean_SF_%i", j), "[0] + [1]*x + [2]*(x*x)",  2, 5 );
		fMean[j]->SetParameters(.1,-.1,-.1);
		fSigma[j] = new TF1( Form("fSigma_SF_%i", j), "[0] + [1]*x + [2]*(x*x)",  2, 5 );
		fSigma[j]->SetParameters(.1,-.1,-.1);
	}

	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );

	// Electron Variables
	electron e;

	// Pion Variables
	std::vector<pion> pi ;
	
	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted
	
	// Set Output file and tree
	TFile * outputFile = new TFile(outFile_name, "RECREATE");
	int nFiles = 0;
	int RunType = 0;
	int inclusive = 0;

	//Set event count array

	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		//create the event reader
		clas12reader c12(files.GetFileName(i).Data());
		auto mcparts = c12.mcparts();		
		int NeventsTotal = c12.getReader().getEntries();       
    	   	int event = 0;	

    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( event%1000000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			evnum  = c12.runconfig()->getEvent();
			runnum = c12.runconfig()->getRun();
			
			///////////////////////////Initialize variables//////////////////////////////////////////////	
			electrons.clear();
			pions.clear();
			pipluses.clear();
			piminuses.clear();
			
			e.Clear();
			pi.clear();

			Ne = Npi = Npips = Npims = 0;

			
			// Get Particles By PID
			electrons   = c12.getByID( 11   );

			Ne      = electrons.size();
			
			if( Ne < 1 ){ continue; } //Keep only events with one electron...
			
			//////////////electron analysis////////////////////
			int e_idx = GetLeadingElectron(electrons, Ne);	
			e.setElectron( Ebeam, electrons[e_idx]);
			////////////////Pion analysis/////////////////
			
				
			if( !anal.applyElectronFiducials( e ) ){ continue; }
			//Electron Min Edep
			if( !anal.applyElectronEDep( e )){continue;}
			if( !anal.applyElectronPCAL( e ) ){continue;}
			
			//Electron SF Cut
			//if( (e.getEpcal() + e.getEecin() + e.getEecout())/e.get3Momentum().Mag() < .17 ){continue;}
			hSF[e.getDC_sector()-1]->Fill( e.get3Momentum().Mag(), (e.getEpcal() + e.getEecin() + e.getEecout())/e.get3Momentum().Mag() );
		
			
		}
	}

	//fill means and sigmas for each 


	std::cout<< "Finished File loop! \n";
	int nBins = hSF[0]->GetXaxis()->GetNbins();
	double bin_centers[6][nBins];
	double means[6][nBins];
	double means_err[6][nBins];
	double sigmas[6][nBins];
	double sigmas_err[6][nBins];
	int filledBins[6] = {0};
	TGraphErrors * gMeans[6];
	TGraphErrors * gSigmas[6];

	for( int sec = 0; sec < 6; sec++ ){
		for( int bin = 1; bin <= nBins; bin++ ){
			if( hSF[sec]->GetXaxis()->GetBinCenter(bin) < 1 ){continue;}
			TH1F * project = (TH1F *)hSF[sec]->ProjectionY("bin", bin, bin);
			if( project->Integral() <= 0 ){ continue; }
			project->Fit("gaus");
			TF1 *gFit = (TF1*)project->GetListOfFunctions()->FindObject("gaus");
			
			means[sec][filledBins[sec]] = gFit->GetParameter(1);
			means_err[sec][filledBins[sec]] = gFit->GetParError(1);
			sigmas[sec][filledBins[sec]] = gFit->GetParameter(2);
			sigmas_err[sec][filledBins[sec]] = gFit->GetParError(2);
			
			bin_centers[sec][filledBins[sec]] = hSF[sec]->GetXaxis()->GetBinCenter(bin);

			filledBins[sec]++;
			
		}
		double xErr[filledBins[sec]] = {0};
		gMeans[sec] = new TGraphErrors( filledBins[sec], bin_centers[sec], means[sec], xErr, means_err[sec] );
		gSigmas[sec] = new TGraphErrors( filledBins[sec], bin_centers[sec], sigmas[sec], xErr, sigmas_err[sec] );
	
		gMeans[sec]->Fit(Form("fMean_SF_%i", sec), "");
		gSigmas[sec]->Fit(Form("fSigma_SF_%i", sec), "");
	}

	std::cout<<"Writing tree to file\n";
   	outputFile->cd();
	for( int j = 0; j < 6; j++ ){
		hSF[j]->Write();
		gMeans[j]->Write(Form("gMean_%i", j));
		gSigmas[j]->Write(Form("gSigma_%i", j));
		fMean[j]->Write();
		fSigma[j]->Write();
	}

	outputFile->Close();

	std::cout<<"Done!\n";

	auto finish = std::chrono::high_resolution_clock::now();
    	std::chrono::duration<double> elapsed = finish - start;

    	std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...................... END OF MAIN.................................//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug ){

	double beam_helicity = p_event->getHelicity();

	bool helFlip = true;
	if      (runnum>=11093 && runnum<=11283)    helFlip = false; // falls, 10.4 GeV period only
	else if (runnum>=11323 && runnum<=11571)    helFlip = false; // winter

	if (helFlip) {
		beam_helicity = -1 * beam_helicity;
	}
	return beam_helicity;
}



int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne){
	// find leading electron as the one with highest energy

	double  leading_e_E = 0;
	int     leading_e_index = 0;


	for (int eIdx=0; eIdx < Ne; eIdx++) {
		double M, P;   

		M = rp[eIdx]->getCalcMass();
		P = rp[eIdx]->par()->getP();
	
		double E_e = sqrt(M*M + P*P);

		if (E_e > leading_e_E) {
			leading_e_index = eIdx;
			leading_e_E     = E_e;
		}
	}

	return leading_e_index;
}




