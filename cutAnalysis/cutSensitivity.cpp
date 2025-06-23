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
//#include "DCfid_SIDIS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
//#include "e_pid.h"
#include "HipoChain.h"
#include "constants.h"
#include "reader.h"

using namespace clas12;
using namespace constants;

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 7 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Beam energy]\n";
		cerr << "	[Run Type] [Single File?] [Inclusive (0, 1)] [Output File Name (no extension)]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
    double Ebeam = atof(argv[2]); // [GeV]
	int RunType = atoi(argv[3]);
	int inclusive =atoi( argv[5]);
	int trial_num =atoi( argv[4]);
	TString outFileName = argv[6]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 

	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	

	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
	analyzer anal(0, torusBending);
	anal.setAnalyzerLevel(RunType);
	anal.loadCutValues(-1, Ebeam);
	anal.loadSamplingFractionParams();
	
	reader runReader;
	runReader.setNumFiles( nFiles);
	runReader.setRunType( RunType );
	runReader.setEnergy( Ebeam );
	
	clas12root::HipoChain files;
    runReader.readRunFiles(files);

	cout<<"Set output files"<<endl;
	
	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );
	double accCharge = 0;
	double beamCurr = 0;
	// Electron Variables
	electron e;

	// Pion Variables
	std::vector<pion> pi ;
	
	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted
	std::vector<int> electronsMC;
	std::vector<std::vector<int>> pionsMC;
	std::vector<std::vector<int>> pidMC;
	
	// Set Output file and tree
	//TFile * outputFile;
	
	//outputFile = new TFile(outFileName + ".root", "RECREATE");
	
	TH1F * ratio_bins[bins_Q2][bins_xB][bins_Z];
	
	for( int q2 = 0; q2 < bins_Q2; q2++ ){
		for( int x = 0; x < bins_xB; x++ ){
			for( int z = 0; z < bins_Z; z++ ){
				ratio_bins[q2][x][z] = new TH1F(Form("ratio_bin_q2_%i_xB_%i_z_%i", q2, x, z), ";R(z);Counts", 250, 0, 1);
			}
		}
	}
	TH1F * helper_1 = new TH1F("helper1", "helper1", bins_Z, Z_min, Z_max);
	TH1F * helper_2 = new TH1F("helper2", "helper2", bins_Z, Z_min, Z_max);
	
	for( int i = 1; i <= bins_Z; i++ ){
		helper_1->SetBinContent(i, -1.);
		helper_2->SetBinContent(i, 4.);
		
		helper_1->SetBinError(i, 0.);
		helper_2->SetBinError(i, 0.);
	}




	for( int tri = 0; tri < 1; tri++ ){
	
		////////////////////////////////////Begin file loop////////////////////////////////////////////////////
		std::cout<<"begin trial : "<<tri<<std::endl;
		TH1F * ratio_pip[bins_Q2][bins_xB];
		TH1F * ratio_pim[bins_Q2][bins_xB];

		for( int q2 = 0; q2 < bins_Q2; q2++ ){
			for( int x = 0; x < bins_xB; x++ ){
				ratio_pip[q2][x] = new TH1F(Form("ratio_pip_q2_%i_xB_%i", q2, x), ";z;Counts", bins_Z, Z_min, Z_max);
				ratio_pim[q2][x] = new TH1F(Form("ratio_pim_q2_%i_xB_%i", q2, x), ";z;Counts", bins_Z, Z_min, Z_max);
			}
		}
		anal.randomizeCuts();
		anal.printCuts();
		
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
				if( RunType > 0 && event%1000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
				if( RunType == 0 && event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
				event++;
				evnum  = c12.runconfig()->getEvent();
				runnum = c12.runconfig()->getRun();
				
				///////////////////////////Initialize variables//////////////////////////////////////////////	
				electrons.clear();
				pions.clear();
				pipluses.clear();
				piminuses.clear();
			
				electronsMC.clear();
				pionsMC.clear();
				
				e.Clear();
				pi.clear();

				Ne = Npi = Npips = Npims = 0;

				/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
				
				// Get Particles By PID
				electrons   = c12.getByID( 11   );
				pipluses    = c12.getByID( 211  );
				piminuses   = c12.getByID(-211  );

				pions = pipluses;
				pions.insert( pions.end(), piminuses.begin(), piminuses.end() );

				Ne      = electrons.size();
				Npi 	= pions.size();
				
				if( Ne < 1 ){ continue; } //Keep only events with one electron...
				if( Npi == 0 && inclusive != 1 ){ continue; }	
			
				//////////////electron analysis////////////////////
				//Find good electrons
				int e_idx = GetLeadingElectron(electrons, Ne);	
				e.setElectron( Ebeam, electrons[e_idx]);
				
				if( !anal.applyElectronDetectorCuts( e )){continue;}
				if( !anal.applyElectronKinematicCuts( e )){continue;}

				////////////////Pion analysis/////////////////
				
				pion pi_dummy;
				//Find good pions			
				for(int i = 0; i < Npi; i++){
					if( inclusive == 1 ){ continue; }
					pi_dummy.Clear();
					pi_dummy.setMCPion( (bool)RunType );
					pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
					if( !anal.applyPionDetectorCuts( pi_dummy, e ) ) {continue;}
					if( !anal.applyPionKinematicCuts( pi_dummy) ) {continue;}


					int charge = pi_dummy.getCharge();
					int this_bin_Q2 = (int)( ( (e.getQ2() - Q2_min)/(Q2_max-Q2_min) )*bins_Q2);
					int this_bin_xB = (int)( ( (e.getXb() - xB_min)/(xB_max-xB_min) )*bins_xB);

					if( charge>0 ){ ratio_pip[this_bin_Q2][this_bin_xB]->Fill( pi_dummy.getZ() );}
					else{ ratio_pim[this_bin_Q2][this_bin_xB]->Fill( pi_dummy.getZ() );}

				}

				
			}

		}

		std::ofstream txtFile_pip, txtFile_pim, txtFile_ratio;
		txtFile_pip.open(outFileName+ "_pip.txt");
		txtFile_pim.open(outFileName+ "_pim.txt");
		txtFile_ratio.open(outFileName+ "_ratio.txt");

		anal.writeCutsToFile(outFileName + "_cuts.txt");
		
		txtFile_pip<< "xB\tQ2\t";
		txtFile_pim<< "xB\tQ2\t";

		for( int z = 1; z <= bins_Z; z++ ){
			txtFile_pip<<z<<"\t";
			txtFile_pim<<z<<"\t";
		}
		txtFile_pip<<"\n";
		txtFile_pim<<"\n";
		for( int j = 1; j <= bins_xB; j++ ){
			for( int i = 1; i <= bins_Q2; i++ ){
				
				txtFile_pip<<j<<"\t"<<i;
				txtFile_pim<<j<<"\t"<<i;

				for( int z = 1; z <= bins_Z; z++ ){
					txtFile_pip<<"\t"<<ratio_pip[i-1][j-1]->GetBinContent(z);
					txtFile_pim<<"\t"<<ratio_pim[i-1][j-1]->GetBinContent(z);
					//if( ratio_pip[i-1][j-1]->GetBinContent(z) < 0 ) continue;
					//ratio_bins[i-1][j-1][z-1]->Fill( ratio_pip[i-1][j-1]->GetBinContent(z) );
				}
				txtFile_pip<<"\n";
				txtFile_pim<<"\n";

			}
		}	
		
		txtFile_pip.close();
		txtFile_pim.close();
		

		std::cout<<"Done!\n";

		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;

		std::cout << "Done. Elapsed time: " << elapsed.count() << std::endl;
	
	}
	
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



