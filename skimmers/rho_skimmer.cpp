/*
 * THINGS TO DO:
 * 1) Write gen-reco matching
 * 2) Figure out Q2, W issue 
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
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "DCfid_SIDIS.h"
#include "SIDISatBAND_auxiliary.h"
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

using namespace cutVals;
using namespace clas12;
SIDISatBAND_auxiliary aux;

const int DC_layers[3] = {6, 18, 36};
const int nRuns[3] = {6, 1, 1};
const int monteCarloRuns[3][6] ={ 
				{7224, 7302, 7304, 7393, 7439, 7520},
				{ 7457 },
				{ 7427 } 
				};

const int nRunsPyth[3] = {1, 1, 1};
const int monteCarloRunsPyth[3][6] ={ 
				{ 7563 },
				{ 7457 },
				{ 7427 } 
				};

void getRunFiles(TString path, TString runList, clas12root::HipoChain &files, int runType, int nFiles, int beamType);

TVector3 RotateVectorTo_qFrame( TVector3 v, TLorentzVector q, TLorentzVector pe );

bool CheckIfElectronPassedSelectionCuts(clashit &e,
					e_pid pid,
                           		TLorentzVector p_e,
                                     	TVector3 v_e,
                                     	int torusBending,
					DCfid_SIDIS dcfid,
					int cuts,
					TH1F ** hists1,
					TH2F ** hists2,
					TH2F ** sector_hists1,
					TH2F ** sector_hists2);

bool CheckIfPionPassedSelectionCuts(int pionCharge,
                                    	clashit pi,
					TLorentzVector p_pi,
					TVector3 v_pi,
					TVector3 v_e,
					int torusBending,
					DCfid_SIDIS dcfid,
					int cuts,
					TH1F ** hists1,
					TH2F ** hists2);
					

TVector3 GetParticleVertex(clas12::region_part_ptr rp);

void SetLorentzVector (TLorentzVector &p4, clas12::region_part_ptr rp);

int GetBeamHelicity( event_ptr p_event, int runnum, int fdebug );

int GetLeadingElectron(std::vector<region_part_ptr> rp, int Ne);

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader * c12);

void ExtractParticleInformation(clashit &e, TLorentzVector &pe, TVector3 &Ve, region_part_ptr rp);

void ExtractParticleInformation(clas12reader * c12, TLorentzVector &par, TVector3 &v_par, int idx);

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList, int job  );

int FindMatch(TLorentzVector p, mcpar_ptr mcparts, std::vector<int> part_list, int PID);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////-------MAIN--------///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv){
			
	auto start = std::chrono::high_resolution_clock::now();

	if( argc < 9 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [# of Files] [Debug mode (0, 1)] [Beam energy]\n";
		cerr << "	[Run Type] [Inclusive (0, 1)] [Output File Name]\n";
		cerr << "	[MC job] [Cut Level (0 - 3)]\n";
		return -1;
	}
	
	int nFiles = atoi(argv[1]); //set 0 to loop over all files,
        int fdebug = atoi(argv[2]);
       	double Ebeam = atof(argv[3]); // [GeV]
	int RunType = atoi(argv[4]);
	int inclusive =atoi( argv[5]);
	TString outFileName = argv[6]; ///volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/10.2/detector_skims/clasdis_7393.root",//, //Enter 
	int job = atoi(argv[7]);
	int cuts = atoi(argv[8]); //0 - all tight, 1 - tight e  && pi vtx, 2 - tight SF, 3 - all loose 
	
	// Check valid beam energy
	if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
		cout<< "Invalid Beam Energy... Set EBeam = 10.2\n"<<endl;
		Ebeam = 10.2;
	}
    	
	// Read cut values
	double torusBending = -1; //outBending = -1, inBending = 1
    	aux.loadCutValues("macros/cuts/BANDcutValues.csv",torusBending);
	auto db = TDatabasePDG::Instance();
	DCfid_SIDIS dcfid;
	e_pid pid;
	pid.setParamsRGB(Ebeam);

	aux.printCutValues();
	
	// Set input file list    
	// Only seems to work with hipo chain method...

	TString DataPath;
	string runList;
	clas12root::HipoChain files;

	setDataPaths( RunType, Ebeam, DataPath, runList, job );	
	getRunFiles( DataPath, runList, files, RunType, nFiles, (int)( ( Ebeam - 10.2 )/0.2 ) );

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
	std::vector<int> pi_sector_DC;//, pi_rich_angle;
	
	std::vector<double> Z, Z_LC, M_x, M_x_d, xF, eta;

	std::vector<region_part_ptr> electrons, pions, pipluses, piminuses; //For reading from hipo file... not outputted

	double M_x_rho, M_rho;	

	// Set output branches
	outTree->Branch("runnum", &runnum);
	outTree->Branch("torus", &torusBending);
	outTree->Branch("evnum", &evnum);
	outTree->Branch("Ebeam", &Ebeam);
	outTree->Branch("beam", &beam);
	outTree->Branch("charge", &charge);
	
	outTree->Branch("Ne", &Ne);
	outTree->Branch("e", &e);
	outTree->Branch("q", &q);
	outTree->Branch("p_e", &p_e);
	outTree->Branch("v_e", &v_e);
		
	outTree->Branch("Npi", &Npi);
	outTree->Branch("Npips", &Npips);
	outTree->Branch("Npims", &Npims);

	outTree->Branch("pi", &pi);
	outTree->Branch("p_pi", &p_pi);
	outTree->Branch("v_pi", &v_pi);
	outTree->Branch("pi_q", &pi_q);
	outTree->Branch("pi_sector_DC", &pi_sector_DC);
	
	outTree->Branch("q2", &Q2);
	outTree->Branch("omega", &omega);
	outTree->Branch("xB", &xB);
	outTree->Branch("y", &y);
	outTree->Branch("w", &W);
	outTree->Branch("W_d", &W_d);

	outTree->Branch("Z", &Z);
	outTree->Branch("Z_LC", &Z_LC);
	outTree->Branch("M_x", &M_x);
	outTree->Branch("M_x_d", &M_x_d);
	outTree->Branch("xF", &xF);
	
	outTree->Branch("M_x_rho", &M_x_rho );
	outTree->Branch("M_rho", &M_rho );

	//Declare skimmer hists
	TH2F ** hSF_no_cuts = new TH2F*[6];
	TH2F ** hSF = new TH2F*[6];

	for( int i = 0; i < 6; i++ ){
		hSF_no_cuts[i] = new TH2F(Form("hSF_no_cuts_%i", i), "hSF_no_cuts", 160, 0, 8, 75, .1, .35);
		hSF[i] = new TH2F(Form("hSF_%i", i), "hSF", 160, 0, 8, 75, .1, .35);
	}

	TH1F * hPositives = new TH1F("hPos", "hPos", 500, 0.75, 1.01);
	TH1F * hNegatives = new TH1F("hMin", "hMin", 500, 0.75, 1.01);
	
	TH2F * hPos_b_p = new TH2F("hPos_b_p", "hPos_b_p", 150, 0, 10, 500, 0.75, 1.01);
	TH2F * hNeg_b_p = new TH2F("hMin_b_p", "hMin_b_p", 150, 0, 10, 500, 0.75, 1.01);

	TString cutlistEl[5] = {"fid", "WV", "minE", "SamplingFraction", "all"};
	TString cutlistPi[4] = {"fid", "vertex", "chi2", "all"};

	TH1F ** hElFid = new TH1F*[5];
	TH1F ** hPiFid = new TH1F*[8];

	TH2F ** hElFid_b_p = new TH2F*[5];
	TH2F ** hPiFid_b_p = new TH2F*[8];

	for( int i = 0; i < 8; i++ ){
		if( i < 5 ){
			hElFid[i] = new TH1F("hEl" + cutlistEl[i], "hEl" + cutlistEl[i], 500, 0.75, 1.01);
			hElFid_b_p[i] = new TH2F("hEl_b_p_" + cutlistEl[i], "hEl_b_p_" + cutlistEl[i], 500, 0, 10, 250, 0.75, 1.01);
		}
		if( i < 4 ){
			hPiFid[i] = new TH1F("hPip" + cutlistPi[i], "hPip" + cutlistPi[i], 500, 0.75, 1.01);
			hPiFid_b_p[i] = new TH2F("hPip_b_p_" + cutlistPi[i], "hPip_b_p_" + cutlistPi[i], 1000, 0, 10, 500, 0.75, 1.01);
		}
		if( i >= 4 ){
			hPiFid[i] = new TH1F("hPim" + cutlistPi[i-4], "hPim" + cutlistPi[i-4], 250, 0.75, 1.01);
			hPiFid_b_p[i] = new TH2F("hPim_b_p_" + cutlistPi[i-4], "hPim_b_p_" + cutlistPi[i-4], 1000, 0, 10, 500, 0.75, 1.01);
		}
	}
	//hElFid[4] = new TH1F( "hVz_loose", "hVz_loose", 1000, -15, 15);


	////////////////////////////////////Begin file loop////////////////////////////////////////////////////
    	for(Int_t i=0;i< files.GetNFiles();i++){//files->GetEntries();i++){
    	    	if (fdebug) std::cout << "reading file " << i+1 <<" of "<<files.GetNFiles()<<"\n ---------------------------------------------------"<< std::endl;
	   	
		//Only skim desired number of files
		if(nFiles != 0 && i > nFiles){break;}	
	
			
		//create the event reader
		//clas12reader c12(files->At(i)->GetTitle(),{0});
		clas12reader c12(files.GetFileName(i).Data());

		int NeventsTotal = c12.getReader().getEntries();       
		cout<<NeventsTotal<<endl;
    	   	int event = 0;
			
    	    	// process the events...
    	    	while((c12.next()==true)){
           		if( RunType > 0 && event%1000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
           		if( RunType == 0 && event%100000 == 0){cout<<"Processing Event: "<<event<< "/"<<NeventsTotal<<endl; }
			event++;
			
			runnum = c12.runconfig()->getRun();
			evnum  = c12.runconfig()->getEvent();
			
			///////////////////////////Initialize variables//////////////////////////////////////////////	
			electrons.clear();
			pions.clear();
			pipluses.clear();
			piminuses.clear();
			charge.clear();
		
			
			e.Clear();
    			p_e = TLorentzVector(0,0,0, db->GetParticle( 11   )->Mass());
    			v_e = TVector3();
    			
			pi.clear();
			p_pi.clear();
			v_pi.clear();

			pi_q.clear();
			M_x_rho = M_rho = 0;	
			
			Ne = Npi = Npips = Npims = 0;
			Q2 = omega = xB = y = W = W_d = 0;			
			pi_sector_DC.clear();

			Z.clear();
			Z_LC.clear();
			M_x.clear();
			M_x_d.clear();
			xF.clear();
			eta.clear();
			
			/////////////////////////////BEGIN EVENT ANALYSIS///////////////////////////
			
			//Fill hists with all negatives and positives
			std::vector<region_part_ptr> allParts = c12.getDetParticles();	
	
			if( RunType == 0 ){
				for( auto particle : allParts ){
					
					if( particle->par()->getCharge() > 0 ){
						hPositives->Fill(particle->getBeta());
						hPos_b_p->Fill(particle->par()->getP(),particle->getBeta() );	
					}
					if( particle->par()->getCharge() < 0 ){
						hNegatives->Fill(particle->getBeta());
						hNeg_b_p->Fill( particle->par()->getP(),particle->getBeta());	
					}

				}
			}

			
				
			// Get Particles By PID
			electrons   = c12.getByID( 11   );
			pipluses    = c12.getByID( 211  );
			piminuses   = c12.getByID(-211  );
			
			pions = pipluses;
			pions.insert( pions.end(), piminuses.begin(), piminuses.end() );
			
			Ne      = electrons.size();
			Npips   = pipluses.size();
			Npims   = piminuses.size();
			Npi 	= pions.size();
			
					
			//if( Ne != 1 ){ continue; } //Keep only events with one electron...
			if( Ne < 1 ){ continue; } //Keep only events with one electron...
			if( Npips < 1 || Npims < 1 ){ continue; }
			if( Npi == 0 && inclusive != 1 ){ continue; }	
			
			//////////////electron analysis////////////////////
			
			//Find good electrons
			int e_idx = GetLeadingElectron(electrons, Ne);	
			ExtractParticleInformation  (e, p_e, v_e, electrons[e_idx]);
			bool selection_e = CheckIfElectronPassedSelectionCuts(e, pid, p_e, v_e, torusBending, dcfid, cuts, hElFid, hElFid_b_p, hSF_no_cuts, hSF);
			e.setSelection(selection_e);				
				
			if(!e.getSelection()){continue;}

			//Compute output quantities

    			q       = beam - p_e;
    			Q2      = -q.Mag2();
			omega   = q.E();
    			xB      = Q2/(2. * aux.Mp * q.E());
			y       = omega / Ebeam;
    			W       = sqrt((p_rest + q).Mag2());
    			W_d     = sqrt((d_rest + q).Mag2());
			

			//Add in electron cuts!


			////////////////Pion analysis/////////////////
			

			//Use dummy variables that will fill vectors
			clashit pi_dummy;
			TLorentzVector p_dummy, pi_q_dummy;
			TVector3 v_dummy, p_pi_q_dummy;
			
			//Use dummy vars for partner in rho decay
			clashit pi_dummy_rho;
			TLorentzVector p_dummy_rho, pi_q_dummy_rho;
			TVector3 v_dummy_rho, p_pi_q_dummy_rho;
			
			bool isGoodRho = false;

			//Find good pions			
			for(int i = 0; i < Npi; i++){
				//Clear dummy variables
				pi_dummy.Clear();
				p_dummy.Clear();
				v_dummy.Clear();
				pi_q_dummy.Clear();
				p_pi_q_dummy.Clear();
				
				pi_dummy_rho.Clear();
				p_dummy_rho.Clear();
				v_dummy_rho.Clear();
				pi_q_dummy_rho.Clear();
				p_pi_q_dummy_rho.Clear();
				
				//Get ith pion info and check if passes cuts
				if( inclusive == 1 ){ continue; }
				ExtractParticleInformation(pi_dummy, p_dummy, v_dummy, pions[i]);

				int charge_temp = (int) pions[i]->par()->getCharge();
				int charge_temp_rho;
				bool selection_pi = CheckIfPionPassedSelectionCuts(charge_temp ,pi_dummy, p_dummy, v_dummy, v_e, torusBending, dcfid, cuts,  hPiFid, hPiFid_b_p); 
				pi_dummy.setSelection(selection_pi);

				if( !selection_pi ){ continue; }
				
				double M_x_temp = (q + p_rest - p_dummy ).Mag();	
				double Z_temp = p_dummy.E()/omega;	


                        	if ( ( M_x_temp < 1.7 || M_x_temp > 5.) ) { continue; }
                        	if ( p_dummy.Vect().Mag() < 1.25 || p_dummy.Vect().Mag() > 5. ) { continue; }
                        	if ( Z_temp < Z_min  ||  Z_temp > Z_max ) { continue; }
                        	if ( p_dummy.Theta()*rad_to_deg < 5. || p_dummy.Theta()*rad_to_deg > 35 ){ continue; }
				

                        	int sector_i = pi_dummy.getDC_sector();
				double theta = p_dummy.Theta()*rad_to_deg;
				double acc_map_pip_min = pips_parameters[sector_i-1][0] + pips_parameters[sector_i-1][1]/p_dummy.Vect().Mag();                      
                        	double acc_map_pim_min = pims_parameters[sector_i-1][0] + pims_parameters[sector_i-1][1]/p_dummy.Vect().Mag();

				if ( theta < acc_map_pip_min && theta < acc_map_pim_min ){ continue; }
		
				//If good pion, search for partner from rho
				for(int j = 0; j < Npi; j++){
					if( isGoodRho ) { break; }; //If good partner found, move on
					charge_temp_rho = (int) pions[j]->par()->getCharge();
					if( charge_temp_rho == charge_temp ){ continue; }
					
					ExtractParticleInformation(pi_dummy_rho, p_dummy_rho, v_dummy_rho, pions[j]);
					bool selection_pi_rho = CheckIfPionPassedSelectionCuts(charge_temp_rho ,pi_dummy_rho, p_dummy_rho, v_dummy_rho, v_e, torusBending, dcfid, cuts,  hPiFid, hPiFid_b_p); 
					if( !selection_pi_rho ){ continue; }
					pi_dummy_rho.setSelection(selection_pi_rho);
					
					M_rho = (p_dummy_rho + p_dummy).Mag();
					M_x_rho = (p_rest + q - p_dummy - p_dummy_rho).Mag();
					
					if( M_x_rho < .8 || M_x_rho > 1.15 ){continue;}
					if( M_rho < .6 || M_rho > .85 ){continue;}
					
					isGoodRho = true;

				}
				

				//If not a good Rho event, keep looking
				if( !isGoodRho ){ continue; }
			
				//If good, fill output variables
			
				pi_q_dummy.SetVectM( RotateVectorTo_qFrame( p_dummy.Vect(), q, p_e ), aux.Mpi );; 					

				Z.push_back(		p_dummy.E()/omega	);
				Z_LC.push_back(		(pi_q_dummy.E() + pi_q_dummy.Pz()) / (q.E() + q.P())	);

				M_x.push_back(		( q + p_rest - p_dummy ).Mag()	);
				M_x_d.push_back(	( q + d_rest - p_dummy ).Mag()	);

				xF.push_back(		2. * (p_dummy.Vect().Dot(q.Vect())) / (q.P() * W)	);

				eta.push_back(		0.5 * log((pi_q_dummy.E()+pi_q_dummy.Pz()) /
								(pi_q_dummy.E()-pi_q_dummy.Pz()))	);

				pi_sector_DC.push_back( pi_dummy.getDC_sector() );	


				pi.push_back(pi_dummy);
				p_pi.push_back(p_dummy);
				v_pi.push_back(v_dummy);
				charge.push_back(charge_temp);
				pi_q.push_back( pi_q_dummy );
				
				//and for partner...
				
				pi_q_dummy_rho.SetVectM( RotateVectorTo_qFrame( p_dummy_rho.Vect(), q, p_e ), aux.Mpi );; 					

				Z.push_back(		p_dummy_rho.E()/omega	);
				Z_LC.push_back(		(pi_q_dummy_rho.E() + pi_q_dummy_rho.Pz()) / (q.E() + q.P())	);

				M_x.push_back(		( q + p_rest - p_dummy_rho ).Mag()	);
				M_x_d.push_back(	( q + d_rest - p_dummy_rho ).Mag()	);

				xF.push_back(		2. * (p_dummy_rho.Vect().Dot(q.Vect())) / (q.P() * W)	);

				eta.push_back(		0.5 * log((pi_q_dummy_rho.E()+pi_q_dummy_rho.Pz()) /
								(pi_q_dummy_rho.E()-pi_q_dummy_rho.Pz()))	);

				pi_sector_DC.push_back( pi_dummy_rho.getDC_sector() );	

				pi.push_back(pi_dummy_rho);
				p_pi.push_back(p_dummy_rho);
				v_pi.push_back(v_dummy_rho);
				charge.push_back(charge_temp_rho);
				pi_q.push_back( pi_q_dummy_rho );


				//If a good set of pions is found, we move to next event	
				break;
			}

			if(p_pi.size() != 2 && inclusive == 0){continue;}
			
	
			outTree->Fill();
			
		}
	
	}
         
   	outputFile->cd();
    	outTree->Write();
	
	hPositives->Write();
	hNegatives->Write();
	hPos_b_p->Write();
	hNeg_b_p->Write();

	//outHistFile->cd();
	for( int i = 0; i < 8; i++ ){
		hPiFid[i]->Write();
		hPiFid_b_p[i]->Write();
	
    	}
	for( int i = 0; i < 5; i++ ){
		hElFid_b_p[i]->Write();
		hElFid[i]->Write();
	
    	}

	
	for( int i = 0; i < 6; i++ ){
		hSF[i]->Write();
		hSF_no_cuts[i]->Write();
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

void setDataPaths( int runType, double eBeam, TString &dataPath, string &runList, int job  ){

	if(runType == 1){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/job_%i/output/%i-",eBeam, job, job);//reco_out/clasdis_rec_", eBeam);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    		//dataPath = "/volatile/clas12/osg/jphelan/job_7393/output/7393-";//reco_out/clasdis_rec_", eBeam);
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/clasdis/%.1f/reco_out/clasdis_rec_", eBeam);
    		//dataPath = "/work/clas12/users/jphelan/SIDIS_at_BAND/generator/dst_";
	}
	
	else if( runType == 2 ){
    		//dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/claspyth/%.1f/proton/reco_out/claspyth_rec_", eBeam);
    		dataPath = "/volatile/clas12/osg/jphelan/job_";
    	}
	
	else if( runType == 3 ){
    		dataPath = Form("/volatile/clas12/users/jphelan/SIDIS/GEMC/claspyth/%.1f/neutron/reco_out/claspyth_rec_", eBeam);
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

	//else{
	//	for( int i = 0; i < 1000; i++){
	//		if( nFiles != 0 && i >= nFiles ) break;
	//		inFile = path + Form("%i.hipo", i+1);
	//		if( gSystem->AccessPathName(inFile) ) continue;
	//		files.Add(inFile.Data());
	//	}
	//}
}


//Check electron detector cuts
bool CheckIfElectronPassedSelectionCuts(clashit &e, 
					e_pid pid, 
					TLorentzVector p_e, 
					TVector3 v_e, 
					int torusBending, 
					DCfid_SIDIS dcfid,
					int cuts,
					TH1F ** hists1,
					TH2F ** hists2,
					TH2F ** sector_hists1,
					TH2F ** sector_hists2){
	// torusBending         torus magnet bending:   ( 1 = inbeding, -1 = outbending    )

	// sometimes the readout-sector is 0
	// Justin B. Estee (June-21): I also had this issue. I am throwing away sector 0. The way you check is plot the (x,y) coordinates of the sector and you will not see any thing. Double check me but I think it is 0.
	hists1[4]->Fill(e.getBeta());
	hists2[4]->Fill( e.getMomentum() ,e.getBeta());
	
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
	
	hists1[0]->Fill(e.getBeta());
	hists2[0]->Fill(e.getMomentum(), e.getBeta());

	if( ! (// fiducial cuts on PCAL
		e.getW() > aux.cutValue_e_PCAL_W
		&&  e.getV() > aux.cutValue_e_PCAL_V)) return false;

	hists1[1]->Fill(e.getBeta());
	hists2[1]->Fill( e.getMomentum() ,e.getBeta());

	if( (cuts == 0 || cuts == 2 || cuts == 3) && ! (
		// Electron Identification Refinement  - PCAL Minimum Energy Deposition Cut
		e.getEpcal() > aux.cutValue_e_E_PCAL)) return false;
		
	hists1[2]->Fill(e.getBeta());
	hists2[2]->Fill(e.getMomentum(),e.getBeta());
	sector_hists1[e.getDC_sector() - 1]->Fill(p_e.P(), (e.getEpcal() + e.getEecin() + e.getEecout())/p_e.P() );

	if( (cuts == 0 || cuts == 1) && !(pid.isElectron(&e)) ) return false;

	if( ! (
		// Sampling fraction cut
		((e.getEpcal() + e.getEecin() + e.getEecout())/p_e.P()) > aux.cutValue_SamplingFraction_min
		&& (e.getEecin()/p_e.P() > aux.cutValue_PCAL_ECIN_SF_min - e.getEpcal()/p_e.P()) // RGA AN puts "<" here mistakenly

		// Cut on z-vertex position: in-bending torus field -13.0 cm < Vz < +12.0 cm
		// Spring 19 and Spring 2020 in-bending.
		// Fall	 2019 (without low-energy-run) was out-bending.
	) ) return false;
	
	if( (cuts == 1 || cuts == 3 ) &&  !((aux.cutValue_Vz_min < v_e.Z()) && (v_e.Z() < aux.cutValue_Vz_max))
	) return false;

	sector_hists2[e.getDC_sector() - 1]->Fill(p_e.P(), (e.getEpcal() + e.getEecin() + e.getEecout())/p_e.P() );


	if( (cuts == 0 || cuts == 2 ) && ! ((v_e.Z() > -5) && (v_e.Z() < -1))
		) return false;
	
	
	hists1[3]->Fill(e.getBeta());
	hists2[3]->Fill( e.getMomentum() ,e.getBeta());
	

	//dcElectrons++;
	return true;
}



// Oo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.
bool CheckIfPionPassedSelectionCuts(int pionCharge, // "pi+" or "pi-"
				clashit pi,
				TLorentzVector p_pi,
				TVector3 v_pi,
				TVector3 v_e,
				int torusBending,
				DCfid_SIDIS dcfid,
				int cuts,
				TH1F ** hists1,
				TH2F ** hists2){

	// decide if pion (pi+ or pi-) passed event selection cuts
	//
	// input:
	// --------
	// DC_x, DC_y   pi drift-chamber coordinates
	// chi2PID      pi chi2PID     (pips_chi2PID)
	// p            pi momentum    (pi.P())
	//
	if( pionCharge > 0 ){
		hists1[3]->Fill(pi.getBeta());
		hists2[3]->Fill(pi.getMomentum(), pi.getBeta());
	}
	else{
		hists1[7]->Fill(pi.getBeta());
		hists2[7]->Fill(pi.getMomentum(),pi.getBeta());
	}

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


	if( pionCharge > 0 ){
		hists1[0]->Fill(pi.getBeta());
		hists2[0]->Fill(pi.getMomentum(), pi.getBeta());
	}
	else{
		hists1[4]->Fill(pi.getBeta());
		hists2[4]->Fill(pi.getMomentum(),pi.getBeta());
	}
	//double chi2PID = pi.getDC_chi2() / pi.getDC_NDF();
	if(! (
	
	// pi+ Identification Refinement - chi2PID vs. momentum
	( aux.Chi2PID_pion_lowerBound( p_pi.P(), C ) < pi.getChi2()
         && pi.getChi2() < aux.Chi2PID_pion_upperBound( p_pi.P() , C ) )
	
       // Cut on the z-Vertex Difference Between Electrons and Hadrons.
	)) { return false; }
	if( pionCharge > 0 ){
		hists1[2]->Fill(pi.getBeta());
		hists2[2]->Fill(pi.getMomentum(), pi.getBeta());
	}
	else{
		hists1[6]->Fill(pi.getBeta());
		hists2[6]->Fill(pi.getMomentum(),pi.getBeta());
	}
       
	if( (cuts == 1 || cuts == 3 ) &&  !( fabs((v_e-v_pi).Z()) < aux.cutValue_Ve_Vpi_dz_max )
	) { return false; }

	if( (cuts == 0 || cuts == 2 ) && !( (v_pi - v_e).Z() > -7 && (v_pi - v_e).Z() < 5 ) ) { return false; }
	
	if( pionCharge > 0 ){
		hists1[1]->Fill(pi.getBeta());
		hists2[1]->Fill(pi.getMomentum(),pi.getBeta());
	}
	else{
		hists1[5]->Fill(pi.getBeta());
		hists2[5]->Fill(pi.getMomentum(),pi.getBeta());
	}

	
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

int GetLeadingElectron(std::vector<int>electrons, int Ne, clas12reader * c12){
	// find leading electron as the one with highest energy
	// generator level
	// note overloaded function
	
	double  leading_e_E = 0;
	int     leading_e_index = 0;


	for (int eIdx=0; eIdx < Ne; eIdx++) {
		double M, P2;   

		//For generator level, need mc particle info
		M = c12->mcparts()->getMass(eIdx);
		P2 = c12->mcparts()->getPx(eIdx)*c12->mcparts()->getPx(eIdx) + 
		c12->mcparts()->getPy(eIdx)*c12->mcparts()->getPy(eIdx) +
		c12->mcparts()->getPz(eIdx)*c12->mcparts()->getPz(eIdx);
	
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
	e.setBeta(rp->par()->getBeta());

	// RICH info
	auto RICH_info = rp->rich();
	double temp_rich_angle(RICH_info->getBest_ch());
	
	if( temp_rich_angle > 0. ){
		e.setRich_beta( 1./(rich_n*cos(temp_rich_angle)) );
	}
	else{ e.setRich_beta(0.); }

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

