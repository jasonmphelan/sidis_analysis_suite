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

using namespace clas12;
using namespace constants;

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

int getThetaBin( double theta );
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

const double binEdges[5] = { 5, 12, 18, 24, 40 };

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 3 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [outFile name (no extension)] [Beam energy (0 for all energies)]\n";
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

	TH1F * hFid_e[2][3][5][6][2]; //[charge][layer][theta bin][sector][weighted?]
	TH1F * hFid_pi[2][3][5][6][2];

	double e_ext[3] = {20, 20,  30};
	double pi_ext[3] = {20, 20,  30};
	

	for( int i = 0; i < 2; i++ ){
	
		for( int j = 0; j < 3; j++ ){
			for( int k = 0; k < 5; k++ ){
				for( int l = 0; l < 6; l++ ){
					hFid_e[i][j][k][l][0] = new TH1F(Form("hFid_e_uw_reg_%i_%i_%i_", j,k,l) + spec[i], ";edge [mm];#chi^2", 100, 0, e_ext[j]);
					hFid_e[i][j][k][l][1] = new TH1F(Form("hFid_e_w_reg_%i_%i_%i_", j,k,l) + spec[i],";edge [mm];#chi^2" , 100, 0, e_ext[j]);
			
					hFid_pi[i][j][k][l][0] = new TH1F(Form("hFid_pi_uw_reg_%i_%i_%i_", j,k,l) + spec[i],";edge [mm];#chi^{2}", 100, 0, pi_ext[j]);
					hFid_pi[i][j][k][l][1] = new TH1F(Form("hFid_pi_w_reg_%i_%i_%i_", j,k,l) + spec[i],";edge [mm];#chi^{2}",  100, 0, pi_ext[j]);
				}
			}
		}
		
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
	TFile * outputFile = new TFile(outFile_name + ".root", "RECREATE");
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
           		if( event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
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

			/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
			//Fill hists with all negatives and positives
			
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
			int e_idx = GetLeadingElectron(electrons, Ne);	
			e.setElectron( Ebeam, electrons[e_idx]);
			////////////////Pion analysis/////////////////
			
			pion pi_dummy;
			//Find good pions			
			for(int i = 0; i < Npi; i++){

				pi_dummy.Clear();
				pi_dummy.setPion( e.getQ(),e.get4Momentum(), pions[i] );
	
			
				int chargeIdx = (int)(pi_dummy.getCharge() < 0);

				//PCAL Fiducials
				if( !anal.applyElectronPCAL( e ) ){continue;}
				
				//Electron Min Edep
				if( !anal.applyElectronEDep( e )){continue;}
				
				//Electron SF Cut
				if( !anal.applyElectronSF( e )){continue;}
				//Electron SF Correcleation
				if( !anal.applyElectronCorrelation( e )){continue;}
				
				int eBin = getThetaBin(e.get3Momentum().Theta()*rad_to_deg);
				int piBin = getThetaBin(pi_dummy.get3Momentum().Theta()*rad_to_deg);
				if( eBin < 0 || piBin < 0 ){continue;}
				if( pi_dummy.getDC_sector() == 0 || e.getDC_sector() == 0 ){continue;}
				if( pi_dummy.getChi2() == 9999 || e.getChi2() == 9999 ){continue;}
				for( int j = 0; j < 3; j++ ){//loop through DC layers
					//cout<<"Pion theta bin : "<<piBin<<std::endl;
					//cout<<"Pion chi2 : "<<pi_dummy.getChi2()<<std::endl;
					//cout<<"Pion edge " << j <<" : "<<pi_dummy.getEdge(j)<<std::endl;
					hFid_e[chargeIdx][j][eBin][e.getDC_sector() - 1][1]->Fill(e.getEdge(j), abs(e.getChi2()));
					hFid_e[chargeIdx][j][eBin][e.getDC_sector() - 1][0]->Fill(e.getEdge(j));	
					hFid_e[chargeIdx][j][0][e.getDC_sector() - 1][1]->Fill(e.getEdge(j), abs(e.getChi2()));
					hFid_e[chargeIdx][j][0][e.getDC_sector() - 1][0]->Fill(e.getEdge(j));
					
					hFid_pi[chargeIdx][j][piBin][pi_dummy.getDC_sector() - 1][1]->Fill(pi_dummy.getEdge(j), abs(pi_dummy.getChi2()));
					hFid_pi[chargeIdx][j][piBin][pi_dummy.getDC_sector() - 1][0]->Fill(pi_dummy.getEdge(j));
					hFid_pi[chargeIdx][j][0][pi_dummy.getDC_sector() - 1][1]->Fill(pi_dummy.getEdge(j), abs(pi_dummy.getChi2()));
					hFid_pi[chargeIdx][j][0][pi_dummy.getDC_sector() - 1][0]->Fill(pi_dummy.getEdge(j));
					
				}	

				//Electron Vertex
				//if( !anal.applyElectronVertex( e )){continue;}
				
				////All Cuts (Check)
				//if( !anal.applyElectronDetectorCuts( e )){continue;}
				//hBeta_p_e[7]->Fill( e.get3Momentum().Mag(), e.getBeta() );
				//hBeta_e[7]->Fill( e.getBeta() );
				//hBeta_p_pi[7][chargeIdx]->Fill( pi_dummy.get3Momentum().Mag(), pi_dummy.getBeta() );
				//hBeta_pi[7][chargeIdx]->Fill( pi_dummy.getBeta() );
				
				//Pion Vertex
				//if( !anal.applyPionDetectorVertex( pi_dummy, e ) ){ continue; }
				

			}

			
		}
	}

	std::cout<< "Finished File loop! \n";
	
	std::cout<<"Writing tree to file\n";
   	outputFile->cd();
	for( int i = 0; i < 2; i++ ){
	
		for( int j = 0; j < 3; j++ ){
			for( int k = 0; k < 5; k++ ){
				for( int l = 0; l < 6; l++ ){
				hFid_e[i][j][k][l][1]->Divide(hFid_e[i][j][k][l][0]);
				hFid_pi[i][j][k][l][1]->Divide(hFid_pi[i][j][k][l][0]);
				
				hFid_e[i][j][k][l][1]->Write();
				hFid_pi[i][j][k][l][1]->Write();
			
				}
			}
		}
		
	}

	outputFile->Close();


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


int getThetaBin( double theta ){
	int bin = -1;

	for( int i = 0; i< 4; i++ ){
		if( theta > binEdges[i] && theta < binEdges[i+1] ){
			return i+1;
		}
	}

	return -1;
}
