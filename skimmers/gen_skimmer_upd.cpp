/*
 * THINGS TO DO:
 * 1) Write gen-reco matching
 * 2) Do comparison to old skimmer
 */

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TChain.h>
//#include <TObject.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "DCfid_SIDIS.h"
//#include "csv_reader.h"
#include "SIDISatBAND_auxiliary.h"
//#include "clashit.cpp"
#include "clashit.h"
#include "e_pid.h"
#include "HipoChain.h"
#include "cut_values.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class vector<TVector3>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class std::vector<clashit>+;
#pragma link C++ class vector<clashit>+;
#endif

using namespace clas12;
using namespace cutVals; 
SIDISatBAND_auxiliary aux;

const int nRuns[3] = {6, 3, 1};
const int monteCarloRuns[3][6] ={ 
				{7224, 7302, 7304, 7393, 7439, 7520},
				{ 7607, 7608, 7769 },
				{ 7427 } 
				};

const int nRunsPyth[3] = {1, 1, 1};
const int monteCarloRunsPyth[3][6] ={ 
				{ 7563 },
				{ 7457 },
				{ 7427 } 
				};

const int DC_layers[3] = {6, 18, 36};

void getRunFiles(TString path, TString runList, clas12root::HipoChain &files, int runType, int nFiles, int beamType);

TVector3 RotateVectorTo_qFrame( TVector3 v, TLorentzVector q, TLorentzVector pe );

bool CheckIfElectronPassedSelectionCuts(clashit &e,
					e_pid pid,
                           		TLorentzVector p_e,
                                     	TVector3 v_e,
                                     	int torusBending,
					DCfid_SIDIS dcfid);

bool CheckIfPionPassedSelectionCuts(int pionCharge,
                                    	clashit pi,
					TLorentzVector p_pi,
					TVector3 v_pi,
					TVector3 v_e,
					int torusBending,
					DCfid_SIDIS dcfid);
					

TVector3 GetParticleVertex(clas12::region_part_ptr rp);

void SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader c12);

void ExtractParticleInformation(clashit &e, TLorentzVector &pe, TVector3 &Ve, region_part_ptr rp);

void ExtractParticleInformation(clas12reader * c12, TLorentzVector &par, TVector3 &v_par, int idx);

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList, int job  );

int FindMatch(TLorentzVector p, mcpar_ptr mcparts, std::vector<int> part_list, int PID);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////


int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 8 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Debug mode (0, 1)] [Beam energy]\n";
		cerr << "	[Run Type] [Inclusive (0, 1)] [Output File Name]\n";
		cerr << "	[MC job] \n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
        int fdebug = atoi(argv[2]);
       	double Ebeam = atof(argv[3]); // [GeV]
	int RunType = atoi(argv[4]);
	int inclusive =atoi( argv[5]);
	TString outFileName = argv[6]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 
	int job = atoi(argv[7]);
			
	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
/*	
	//Check single file
	if( singleFile != 0 ||singleFile != 1 ){
		cout<<"Invalid 'single file' parameter.... break"<<endl;
		break;	
	}
*/
    	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
    	aux.loadCutValues(torusBending);
	auto db = TDatabasePDG::Instance();
	DCfid_SIDIS dcfid;
	e_pid pid;
	pid.setParamsRGB(Ebeam);
	
	// Set input file list    
	// Only seems to work with hipo chain method...
	
	TString DataPath;
	string runList;
	clas12root::HipoChain files;

	setDataPaths( RunType, Ebeam, DataPath, runList, job );	
	getRunFiles( DataPath, runList, files, RunType, nFiles, (int)( (Ebeam - 10.2)/.2 ) );

	// Set Output file and tree
	TFile * outputFile = new TFile(outFileName, "RECREATE");
	TTree * outTree = new TTree("ePi", "(e,e'pi) event  information");
	
	// Set output variables
	int Ne,Npi, Npips, Npims, runnum, evnum;
	double torus_setting;
	TLorentzVector beam( 0, 0, Ebeam, Ebeam );

	TLorentzVector p_rest( 0, 0, 0, aux.Mp );
	TLorentzVector d_rest( 0, 0, 0, aux.Md );

	// Electron Variables
	double Q2, omega, xB, y, W, W_d;
	TLorentzVector q;			
	clashit e;
	TLorentzVector p_e;
	TVector3 v_e;
	
	// Pion Variables
	std::vector<clashit> pi ;
	std::vector<TLorentzVector> p_pi, pi_q;
	std::vector<TVector3> v_pi;
	std::vector<int> charge;
	std::vector<int> pi_sector_DC;
	
	std::vector<double> Z, Z_LC, M_x, M_x_d, xF, eta;

	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted
	
	//Need to save relevant Generator information for GEMC sim
	int NeGen, NpiGen, NpipsGen, NpimsGen;
	double Q2_gen, omega_gen, xB_gen, y_gen, W_gen;
	std::vector<double> Z_gen, M_x_gen;	
	TLorentzVector p_e_gen, q_gen;
	std::vector<TLorentzVector> p_pi_gen, p_parent;

	std::vector<int> electronsMC, pionsMC, parentPID;

	double weight;

	// Set output branches
	outTree->Branch("runnum", &runnum);
	outTree->Branch("torus", &torusBending);
	outTree->Branch("evnum", &evnum);
	outTree->Branch("Ebeam", &Ebeam);
	outTree->Branch("charge", &charge);
	
	outTree->Branch("Ne", &Ne);
	//outTree->Branch("e", &e);
	//outTree->Branch("p_e", &p_e);
	//outTree->Branch("v_e", &v_e);
		
	outTree->Branch("Npi", &Npi);
	outTree->Branch("Npips", &Npips);
	outTree->Branch("Npims", &Npims);

	//outTree->Branch("pi", &pi);
	//outTree->Branch("p_pi", &p_pi);
	//outTree->Branch("v_pi", &v_pi);
	//outTree->Branch("pi_q", &pi_q);
	outTree->Branch("pi_sector_DC", &pi_sector_DC);
	
	outTree->Branch("Ne", &NeGen);
	outTree->Branch("p_e", &p_e_gen);
	
	outTree->Branch("Npi", &NpiGen);
	outTree->Branch("Npips", &NpipsGen);
	outTree->Branch("Npims", &NpimsGen);

	outTree->Branch("p_pi", &p_pi_gen);

	//outTree->Branch("Q2", &Q2);
	//outTree->Branch("omega", &omega);
	//outTree->Branch("xB", &xB);
	//outTree->Branch("y", &y);
	//outTree->Branch("W", &W);
	//outTree->Branch("W_d", &W_d);

	//outTree->Branch("Z", &Z);
	//outTree->Branch("Z_LC", &Z_LC);
	//outTree->Branch("M_x", &M_x);
	//outTree->Branch("M_x_d", &M_x_d);
	//outTree->Branch("xF", &xF);
	
	outTree->Branch("q2", &Q2_gen);
	outTree->Branch("omega", &omega_gen);
	outTree->Branch("xB", &xB_gen);
	outTree->Branch("y", &y_gen);
	outTree->Branch("w", &W_gen);
	
	outTree->Branch("Z", &Z_gen);
	outTree->Branch("M_x", &M_x_gen);

	outTree->Branch("parentPID", &parentPID);
	outTree->Branch("p_parent", &p_parent);

	outTree->Branch("weight", &weight);

	int goodElectron = 0;

	//Begin file loop
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	if (fdebug) std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;
	   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
				
		//create the event reader
		//clas12reader c12(files->At(i)->GetTitle(),{0});
		clas12reader c12(files.GetFileName(i).Data());
		auto mcparts = c12.mcparts();		

		int NeventsTotal = c12.getReader().getEntries();       
		cout<<NeventsTotal<<endl;
    	   	int event = 0;
			
    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( RunType > 0 && event%10000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
           		if( RunType == 0 && event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			
			//Get run and event info	
			runnum = c12.runconfig()->getRun();
			evnum  = c12.runconfig()->getEvent();
			
			//Initialize variables	
			electrons.clear();
			pions.clear();
			pipluses.clear();
			piminuses.clear();
			charge.clear();
		
			
			e.Clear();
    			p_e = TLorentzVector(0,0,0, db->GetParticle( 11   )->Mass());
    			v_e = TVector3();
    			
			p_e_gen = TLorentzVector(0,0,0, db->GetParticle( 11   )->Mass());
			q_gen = TLorentzVector( 0, 0, 0, 0 );

			pi.clear();
			p_pi.clear();
			v_pi.clear();

			pi_q.clear();
		
			p_pi_gen.clear();
			
			Ne = Npi = Npips = Npims = 0;
			NeGen = NpiGen = NpipsGen = NpimsGen = 0;
			Q2 = omega = xB = y = W = W_d = 0;			
			Q2_gen = omega_gen = xB_gen = y_gen = W_gen = 0;			
			pi_sector_DC.clear();

			electronsMC.clear();
			pionsMC.clear();

			Z.clear();
			Z_LC.clear();
			M_x.clear();
			M_x_d.clear();
			xF.clear();
			eta.clear();
			
			Z_gen.clear();
			M_x_gen.clear();

			p_parent.clear();
			parentPID.clear();

			weight = 0;
			/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
			// Get Particles By Type
		
			weight = c12.mcevent()->getWeight();
			int nMcPart = c12.mcevent()->getNpart();
			
			int mcId;
			for ( int i = 1; i<nMcPart; i++ ){
				mcId = c12.mcparts()->getPid(i);
				switch (mcId){
					case 11:
						electronsMC.push_back(i);
						break;
					case 211:
						Npips++;
						pionsMC.push_back(i);
						break;
					case -211:
						Npims++;
						pionsMC.push_back(i);
						break;
				}
			}
		
			Ne = electronsMC.size();
			Npi = pionsMC.size();
			if( Ne < 1 ){ continue; } //Keep only events with one electron...
			if( Npi == 0 && inclusive != 1 ){ 
			//	cout<<"NO GOOD PIONS"<<endl;
				continue; }	
			if( Npi != Npips + Npims ){ 
				cout<<"ERROR IN COMBINING PIPLUS AND PIMINUS!!!\n";
				break;
			}
			goodElectron++;
			
			//////////////electron analysis////////////////////
		
	
			//Find and fill good electrons

			int e_idx = electronsMC[0];
			//int e_idx = GetLeadingElectron(electronsMC, Ne, c12);	
			//ExtractParticleInformation  (e, p_e, v_e, electrons[e_idx]);
			//bool selection_e = CheckIfElectronPassedSelectionCuts(e, pid, p_e, v_e, torusBending, dcfid);
			//e.setSelection(selection_e);				
	
			//if(!e.getSelection()){continue;}
			
			//Compute output quantities

    			//q       = beam - p_e;
    			//Q2      = -q.Mag2();
    			//omega   = q.E();
    			//xB      = Q2/(2. * aux.Mp * q.E());
    			//y       = omega / Ebeam;
    			//W       = sqrt((p_rest + q).Mag2());
    			//W_d     = sqrt((d_rest + q).Mag2());
			
			//Fill generator values (note, match for electron unambiguous)
			//If no good electron, save time by not doing this
			mcparts->setEntry(e_idx);
								
			p_e_gen.SetPxPyPzE(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), sqrt( pow( mcparts->getP(), 2 ) + pow( mcparts->getMass(), 2 ) ) );
			q_gen = beam - p_e_gen;
				
			//q_gen.SetPxPyPzE(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), sqrt( pow( mcparts->getP(), 2 ) + pow( mcparts->getMass(), 2 ) ) );
			Q2_gen      = -q_gen.Mag2();
			omega_gen   = q_gen.E();
    			xB_gen      = Q2_gen/(2. * aux.Mp * q_gen.E());
    			y_gen       = omega_gen / Ebeam;
    			W_gen       = sqrt((p_rest + q_gen).Mag2());
			
			
			//if( W_gen < 2.5 ) { continue; }
			//if( Q2_gen < Q2_min || Q2_gen > Q2_max ) { continue; }
                	//if( xB_gen < xB_min || xB_gen > xB_max ) { continue; }
                	//if( y_gen > 0.75 ) { continue; }
                	//if( p_e_gen.Vect().Mag() <3. || p_e_gen.Vect().Mag() > 10.6 ) { continue; }
                	//if( p_e_gen.Theta()*rad_to_deg < 5. || p_e_gen.Theta()*rad_to_deg > 35. ){ continue; }

			////////////////Pion analysis/////////////////
			
			//Find and fill good pions
	

			//Use dummy variables that will fill vectors
			clashit pi_dummy;
			TLorentzVector p_dummy, pi_q_dummy;
			TVector3 v_dummy, p_pi_q_dummy;
			
	
			for(int i = 0; i < Npi; i++){
				if( inclusive == 1 ){ continue; }
				
				mcparts->setEntry(pionsMC[i]);
				
				TLorentzVector pi_gen_dummy(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), sqrt( pow( mcparts->getP(), 2 ) + pow( mcparts->getMass(), 2 ) ) );

				double M_x_temp = (q_gen + p_rest - pi_gen_dummy).Mag();                       
				double Z_temp = pi_gen_dummy.E()/omega_gen; 

				//do cuts
				
				//if ( M_x_temp < 1.7 || M_x_temp > 5. ) { continue; }
                        	//if ( M_x[i] < 1.5 || M_x[i] > 5. ) { continue; }
                        	//if ( pi_gen_dummy.Vect().Mag() < 1.25 || pi_gen_dummy.Vect().Mag() > 5. ) { continue; }
                        	//if ( Z_temp < Z_min  ||  Z_temp > Z_max ) { continue; }
                        	//if ( pi_gen_dummy.Theta()*rad_to_deg < 5. || pi_gen_dummy.Theta()*rad_to_deg > 35 ){ continue; }
				

				double phi = pi_gen_dummy.Vect().Phi();

				int charge_temp = (int) ( (double) mcparts->getPid() / 211. );  
				int sector_i = 0;

				if( charge_temp > 0 ){    
                                        if( phi > -0.8 && phi < 0.25 ){ sector_i = 1; }
                                        else if( phi >= 0.25 && phi < 1.3 ){ sector_i = 2; }
                                        else if( phi >= 1.3 && phi <= 2.35 ){ sector_i = 3; }
                                        else if( phi > 2.35 || phi < -2.9  ){ sector_i = 4; }
                                        else if( phi > -2.9 && phi < -1.85){ sector_i = 5; }
                                        else{ sector_i = 6; }
                                }
                                if( charge_temp < 0 ){
                                        if( phi > -0.25 && phi < 0.8 ){ sector_i = 1; }
                                        else if( phi >= 0.8 && phi < 1.85 ){ sector_i = 2; }
                                        else if( phi >= 1.85 && phi <= 2.9 ){ sector_i = 3; }
                                        else if( phi > 2.9 || phi < -2.4  ){ sector_i = 4; }
                                        else if( phi > -2.4 && phi < -1.25){ sector_i = 5; }
                                        else{ sector_i = 6; }
                                }
			
				if( sector_i == 0 ) continue;

				pi_sector_DC.push_back( sector_i );
	
				//fill vectors
				
				charge.push_back( charge_temp );
				
				p_pi_gen.push_back( pi_gen_dummy );
				Z_gen.push_back( pi_gen_dummy.E()/omega_gen );
				M_x_gen.push_back( (q_gen + p_rest - pi_gen_dummy).Mag() );
				
				int parentIdx = mcparts->getParent() - 1;
					
				mcparts->setEntry( parentIdx );
				parentPID.push_back( mcparts->getPid() );
					
				TLorentzVector parent_dummy(mcparts->getPx(), mcparts->getPy(), mcparts->getPz(), sqrt( pow( mcparts->getP(), 2 ) + pow( mcparts->getMass(), 2 ) ) );
				 	
				p_parent.push_back( parent_dummy );

				
				//Clear dummy variables
				pi_dummy.Clear();
				p_dummy.Clear();
				v_dummy.Clear();
				pi_q_dummy.Clear();
				p_pi_q_dummy.Clear();
			}

			if(p_pi_gen.size() == 0 && inclusive == 0){continue;}
			
			outTree->Fill();
			
		}
	
	}
         	
 	cout<<"Number of electrons : "<<goodElectron<<endl;
	
   	outputFile->cd();
    	outTree->Write();
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

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList, int job  ){

	if(runType == 1){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/reco_out/clasdis_rec_", eBeam);
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/job_%i/output/%i-",eBeam, job, job);//reco_out/clasdis_rec_", eBeam);
    		//dataPath = Form("/volatile/clas12/osg/jphelan/job_%i/output/%i-", job, job);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    		//dataPath = "/work/clas12/users/jphelan/SIDIS_at_BAND/generator/dst_";
    		//dataPath = "/work
	}
	
	else if( runType == 2 ){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/claspyth/%.1f/proton/job_7177/output/7177-", eBeam);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    	}
	
	else if( runType == 3 ){
    		dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/claspyth/%.1f/neutron/job_7216/output/7216-", eBeam);
    	}
	

	else{ //runType == 0 is data
		TString path_temp;
	
		if(eBeam == 10.6){ 
			runList = "macros/runlists/good_runs_10-6.txt";
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else if(eBeam == 10.4){ 
			runList = "macros/runlists/good_runs_10-4.txt";
			path_temp = "spring2020/torus-1/pass2/v1/dst/train/sidisdvcs/sidisdvcs_";
		}
	
		else{
			runList = "macros/runlists/good_runs_10-2-final.txt";
			path_temp = "spring2019/torus-1/pass2/v0/dst/train/sidisdvcs/sidisdvcs_";
    		}
    		
		dataPath = "/cache/clas12/rg-b/production/recon/"+path_temp;	
    	}
}

void getRunFiles(TString path, TString runList, clas12root::HipoChain &files, int runType, int nFiles, int beamType){
	std::ifstream stream;
	stream.open(runList);
	string runNum;
	TString inFile;

	if(runType == 0){
		int i = 0;
		while(std::getline(stream, runNum)){
			if(nFiles != 0 && i >= nFiles ) break;
			TString run(runNum);
			inFile = path + run+".hipo";
			//inFile = path + "006421.hipo";
			
			cout<<inFile<<endl;
			files.Add(inFile.Data());		
			i++;
		}
	}
	else if (runType == 1){
		for( int j = 0; j < nRuns[beamType]; j++ ){
			for( int i = 0; i < 75000; i++){
				if( nFiles != 0 && i >= nFiles ) break;
				inFile = path + Form("%i/output/%i-%i.hipo", monteCarloRuns[beamType][j], monteCarloRuns[beamType][j], i+1);
				if( gSystem->AccessPathName(inFile) ) continue;
				files.Add(inFile.Data());
			}
		}
	}
	else {
		for( int j = 0; j < nRunsPyth[beamType]; j++ ){
			for( int i = 0; i < 75000; i++){
				if( nFiles != 0 && i >= nFiles ) break;
				inFile = path + Form("%i/output/%i-%i.hipo", monteCarloRunsPyth[beamType][j], monteCarloRunsPyth[beamType][j], i+1);
				if( gSystem->AccessPathName(inFile) ) continue;
				files.Add(inFile.Data());
			}
		}
	}
}


//Check electron detector cuts
bool CheckIfElectronPassedSelectionCuts(clashit &e, 
					e_pid pid, 
					TLorentzVector p_e, 
					TVector3 v_e, 
					int torusBending, 
					DCfid_SIDIS dcfid){
	// torusBending         torus magnet bending:   ( 1 = inbeding, -1 = outbending    )

	// sometimes the readout-sector is 0
	// Justin B. Estee (June-21): I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing. Double check me but I think it is 0.
	if (e.getDC_sector() == 0) return false;

	double e_DC_x[3] = {e.getDC_x1(), e.getDC_x2(), e.getDC_x3()};
	double e_DC_y[3] = {e.getDC_y1(), e.getDC_y2(), e.getDC_y3()};
	double e_DC_z[3] = {e.getDC_z1(), e.getDC_z2(), e.getDC_z3()};



	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)

		int bending  = 1 ? (torusBending==-1) : 0;
		bool DC_fid  = dcfid.DC_fid_xy_sidis( 11,                 // particle PID,
						e_DC_x[regionIdx],  // x
						e_DC_y[regionIdx],  // y
						e.getDC_sector(),        // sector
						regionIdx+1,        // layer
						bending );           // torus bending
		if (DC_fid == false) { return false; }
	}	


	if( ! (// fiducial cuts on PCAL
		e.getW() > aux.cutValue_e_PCAL_W
		&&  e.getV() > aux.cutValue_e_PCAL_V

		// Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
		&&  e.getEpcal() > aux.cutValue_e_E_PCAL

		// Sampling fraction cut
		&& ((e.getEpcal() + e.getEecin() + e.getEecout())/p_e.P()) > aux.cutValue_SamplingFraction_min
		&& (e.getEecin()/p_e.P() > aux.cutValue_PCAL_ECIN_SF_min - e.getEpcal()/p_e.P()) // RGA AN puts "<" here mistakenly

		// Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
		// Spring 19 and Spring 2020 in-bending.
		// Fall	 2019 (without low-energy-run) was out-bending.
		//&&  ((aux.cutValue_Vz_min < v_e.Z()) && (v_e.Z() < aux.cutValue_Vz_max))

		&&  ((v_e.Z() > -5) && (v_e.Z() < -1))
		)) return false;
	
	if( !(pid.isElectron(&e)) ) return false;

	return true;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(int pionCharge, // "pi+" or "pi-"
				clashit pi,
				TLorentzVector p_pi,
				TVector3 v_pi,
				TVector3 v_e,
				int torusBending,
				DCfid_SIDIS dcfid){

	// decide if pion (pi+ or pi-) passed event selection cuts
	//
	// input:
	// --------
	// DC_x, DC_y   pi drift-chamber coordinates
	// chi2PID      pi chi2PID     (pips_chi2PID)
	// p            pi momentum    (pi.P())
	//

	if (pi.getDC_sector() == 0) { return false;}
	
	int PDGcode;
	double    C;
	
	if (pionCharge > 0){
		PDGcode = 211;
		C       = 0.88;
	} 
	else if (pionCharge < 0) {
		PDGcode = -211;
		C       = 0.93;
	} 
	else {
		std::cout << "π charge ill-defined, returning false" << std::endl;
		return false;
	}

	double DC_x[3] = {pi.getDC_x1(), pi.getDC_x2(), pi.getDC_x3()};
	double DC_y[3] = {pi.getDC_y1(), pi.getDC_y2(), pi.getDC_y3()};
	double DC_z[3] = {pi.getDC_z1(), pi.getDC_z2(), pi.getDC_z3()};

	for (int regionIdx=0; regionIdx<3; regionIdx++) {
		// DC_e_fid:
		// sector:  1-6
		// layer:   1-3
		// bending: 0(out)/1(in)
		
		int bending  = 1 ? (torusBending==-1) : 0;
		bool DC_fid  = dcfid.DC_fid_th_ph_sidis(PDGcode,            // particle PID
							DC_x[regionIdx],    // x
							DC_y[regionIdx],    // y
							DC_z[regionIdx],    // z
							pi.getDC_sector(),          // sector
							regionIdx+1,        // layer
							bending);           // torus bending
		
		if (DC_fid == false) { return false; }
	}


	//double chi2PID = pi.getDC_chi2() / pi.getDC_NDF();
	if(! (
	
	// pi+ Identification Refinement - chi2PID vs. momentum
	( aux.Chi2PID_pion_lowerBound( p_pi.P(), C ) < pi.getChi2()
         && pi.getChi2() < aux.Chi2PID_pion_upperBound( p_pi.P() , C ) )
       
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
       
	&&  ( fabs((v_e-v_pi).Z()) < aux.cutValue_Ve_Vpi_dz_max ))
	) { return false; }
	
	return true;
}



TVector3 GetParticleVertex(clas12::region_part_ptr rp){
	TVector3 V(rp->par()->getVx(),
			rp->par()->getVy(),
			rp->par()->getVz());
	return V;
}

void SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp){

	p4.SetXYZM(rp->par()->getPx(),
			rp->par()->getPy(),
			rp->par()->getPz(),
			p4.M());
}



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

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader c12){
	// find leading electron as the one with highest energy
	// generator level
	// note overloaded function
	
	double  leading_e_E = 0;
	int     leading_e_index = 0;


	for (int eIdx=1; eIdx < Ne; eIdx++) {
		double M, P2;   

		//For generator level, need mc particle info
		M = c12.mcparts()->getMass(eIdx);
		P2 = c12.mcparts()->getPx(eIdx)*c12.mcparts()->getPx(eIdx) + 
		c12.mcparts()->getPy(eIdx)*c12.mcparts()->getPy(eIdx) +
		c12.mcparts()->getPz(eIdx)*c12.mcparts()->getPz(eIdx);
	
		double E_e = sqrt(M*M + P2);

		if (E_e > leading_e_E) {
			leading_e_index = eIdx;
			leading_e_E     = E_e;
		}
	}
	return leading_e_index;
}


void ExtractParticleInformation(clashit &e, TLorentzVector &pe, TVector3 &Ve, region_part_ptr rp){
	// set leading electron 4-momentum
	SetLorentzVector(pe , rp);
	// set leading electron vertex
	Ve              = GetParticleVertex( rp );
	
	e.setMomentum( pe.Vect().Mag() );
	
	e.setCharge( rp->par()->getCharge() );

	e.setChi2(rp->par()->getChi2Pid());

	// detector information on electron
	auto e_PCAL_info= rp->cal(PCAL);
	e.setEpcal(e_PCAL_info->getEnergy());
	e.setSector(e_PCAL_info->getSector());
	e.setV( e_PCAL_info->getLv());
	e.setW(e_PCAL_info->getLw());
	e.setEecin(rp->cal(ECIN)->getEnergy());
	e.setEecout(rp->cal(ECOUT)->getEnergy());

	// hit position in PCAL
	e.setPCal_X(e_PCAL_info->getX());
	e.setPCal_Y( e_PCAL_info->getY());
	e.setPCal_Z(e_PCAL_info->getZ());

	// Sampling Fraction		
	e.setEoP((e.getEpcal() + e.getEecin() + e.getEecout())/pe.P());

	// Drift Chamber tracking system
	auto e_DC_info  = rp->trk(DC);
	e.setDC_sector(e_DC_info->getSector()); // tracking sector
	e.setDC_chi2(e_DC_info->getChi2());  // tracking chi^2/NDF
	e.setDC_NDF(e_DC_info->getNDF());  // tracking chi^2/NDF

	e.setDC_x1(rp->traj(DC,DC_layers[0])->getX());
	e.setDC_y1(rp->traj(DC,DC_layers[0])->getY());	
	e.setDC_z1(rp->traj(DC,DC_layers[0])->getZ());

	e.setDC_x2(rp->traj(DC,DC_layers[1])->getX());
	e.setDC_y2(rp->traj(DC,DC_layers[1])->getY());	
	e.setDC_z2(rp->traj(DC,DC_layers[1])->getZ());

	e.setDC_x3(rp->traj(DC,DC_layers[2])->getX());
	e.setDC_y3(rp->traj(DC,DC_layers[2])->getY());	
	e.setDC_z3(rp->traj(DC,DC_layers[2])->getZ());

}

//Generator Level
void ExtractParticleInformation(clas12reader * c12, TLorentzVector &par, TVector3 &v_par, int idx){
	par.SetXYZM(c12->mcparts()->getPx(idx),
		c12->mcparts()->getPy(idx),
		c12->mcparts()->getPz(idx),
		c12->mcparts()->getMass(idx));

	v_par.SetXYZ(c12->mcparts()->getVx(idx),
		c12->mcparts()->getVy(idx),
		c12->mcparts()->getVz(idx));
}

TVector3 RotateVectorTo_qFrame( TVector3 v, TLorentzVector q, TLorentzVector pe ){
    	// move to q-Pe system: q is the z axis, Pe is in x-z plane: Pe=(Pe[x],0,Pe[q])
    	v.RotateZ( -q.Phi()  );
    	v.RotateY( -q.Theta() );
    	v.RotateZ( -pe.Phi() );

	return v;
}


int FindMatch(TLorentzVector p, mcpar_ptr mcparts, std::vector<int> part_list, int PID){

	double temp_min_dp = 9999;
	double p_idx = -1;
	double temp_min_dtheta = 9999;
	double theta_idx = -1;

	for ( int i : part_list ){
		mcparts->setEntry(i);

		if( mcparts->getPid() != PID ){ continue; }		

		if( abs( p.Vect().Mag() - mcparts->getP() ) < temp_min_dp ){	
			temp_min_dp = abs( p.Vect().Mag() - mcparts->getP() );
			p_idx = i;
		}		
		if( abs( p.Theta() - mcparts->getTheta() ) < temp_min_dtheta ){	
			temp_min_dtheta = abs( p.Theta() - mcparts->getTheta() );
			theta_idx = i;
		}	
	}
	
	if( theta_idx != p_idx ){
		return -1;
	}
	
	return p_idx;
}

