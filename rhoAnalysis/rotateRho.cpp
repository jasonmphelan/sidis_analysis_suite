#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <TROOT.h>
#include <chrono>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TEventList.h"
#include "TRandom3.h"
#include "CLAS.h"
#include "electron.h"
#include "pion.h"
#include "analyzer.h"
#include "constants.h"
#include <thread>
using namespace constants;

#define CORR_PATH _DATA

using std::cerr;
using std::isfinite;
using std::cout;
using std::endl;
using std::ofstream;
using std::thread;

const int thread_MAX = 50;

struct eventData{
	TLorentzVector e4;
	TLorentzVector q;
	pion pi_1;
	pion pi_2;
	bool isGoodPion_1;
	bool isGoodPion_2;
};
struct analData{
	analyzer * anal;
	TRandom3 * gen;
	double err_level;
	int acc_match;
	TString in_name;
	TString out_name;
	int chunkSize;
};

std::array<double, 6> rotateEvent(eventData this_event, analData this_anal);
TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TVector3 pi_q );
void fill_chunk( TTree * inTree, TTree** tree_chunk, int chunkSize, int i_thread);
void rotate_chunk( analData this_anal, int i_thread);


int main( int argc, char** argv){

	ROOT::EnableThreadSafety();

	if( argc < 7 ){
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./code [Input File] [Output File] [beam energy] [Rel Uncertainty Level (%)] \n";
		return -1;
	}

	TString in_name = argv[1];
    TString out_name = argv[2];
	double energy = atof( argv[3] );
	double err_level = atof( argv[4] );
	int acc_match = atoi( argv[5] );
	TString acc_name = argv[6];
	int i_thread = atoi( argv[7] );

	//Set analysis basics

	TRandom3 * gen = new TRandom3(0);

	analyzer * anal = new analyzer(0, -1);
	anal->loadAcceptanceMap( (TString)_DATA + (TString)"/acceptance_map/"+acc_name);//%.1f.root", energy));
	anal->loadMatchingFunctions();

// Number of threads to use
int numThreads = thread_MAX;//thread::hardware_concurrency();
int chunkSize  = 500;//event_total / (numThreads-1);
cout<<"Chunk Size : "<<chunkSize<<std::endl;

	analData this_anal = {anal, gen, err_level, acc_match, in_name, out_name, chunkSize};

	//Open input file/tree
    TFile * file_rec = new TFile(in_name);
	TTree * inTree = (TTree*)file_rec->Get("ePi");
	int event_total = inTree->GetEntries();
	
	

	//Create chunky trees
	//TTree ** tree_chunk = new TTree*[numThreads];
	//TTree ** out_chunk =new TTree*[numThreads];
	//TFile * outFile = new TFile(out_name + ".root", "RECREATE");
	//for( int chunk=i_thread; chunk < i_thread+1; chunk++ ){
	//	tree_chunk[chunk] = inTree->CloneTree(0);
	//	out_chunk[chunk] =  new TTree("ePi", "(e,e'pi) event  information");
	//}

	//Set Branch addresses
	/*
	electron * e = nullptr;
	TLorentzVector * beam = nullptr;
	std::vector<pion>* pi = nullptr;
	std::vector<bool>* isGoodPion = nullptr;
	std::vector<bool>* isGoodPion_acc = nullptr;
	double Mx_2pi;
	double M_rho;

	inTree->SetBranchAddress("e", &e);
	inTree->SetBranchAddress("Mx_2pi", &Mx_2pi);
	inTree->SetBranchAddress("M_rho", &M_rho);
	inTree->SetBranchAddress("pi", &pi);
	inTree->SetBranchAddress("isGoodPion_no_acc", &isGoodPion);
	inTree->SetBranchAddress("isGoodPion", &isGoodPion_acc);
	inTree->SetBranchAddress("beam", &beam);
	*/
	//Fill chunks
	//fill_chunk(inTree, tree_chunk, chunkSize, i_thread);

	//cout<<"Split trees complete\n";
	
	std::vector<std::thread> threads;
	for( int i = 27; i < 28; i++ ){

		threads.emplace_back(rotate_chunk, this_anal, i);
		
	}
	for (auto& t : threads) t.join();
	//TFile f(Form(out_name + "_chunk_%i.root", i_thread), "RECREATE");
	//out_chunk[i_thread]->Write();
	//f.Close();
	//file_rec->Close();
}


std::array<double, 6> rotateEvent(eventData this_event, analData this_anal ){
	
	bool detectedPion[2] = {true, true};
	bool goodPion[2] = {true, true};
	bool accPion[2] = {true, true};

	double onePiEvents[2] = {0};
	double twoPiEvents = 0;
	double trials = 0;

	double rhoWeight[2];
	double corr_err[2];
	rhoWeight[0] = 0; // (pi+, pi-)
	rhoWeight[1] = 0; // (pi+, pi-)
	corr_err[0] = 999;
	corr_err[1] = 999;

	//Keep rotating while (err > input) for each pion.  Also set a mininum number of trials
	while( (this_event.isGoodPion_1 && corr_err[0] > this_anal.err_level/100.) || (this_event.isGoodPion_2 && corr_err[1] > this_anal.err_level/100.) || trials < 1000){
		trials++;
		if( trials > 500000 ) break;

		//if( trials > 50000 && ( corr_err[0] > 900 || corr_err[1] > 900 ) ){ break; }
		detectedPion[0] = false;
		detectedPion[1] = false;
		goodPion[0] = false;
		goodPion[1] = false;
		accPion[0] = false;
		accPion[1] = false;

		//initialize rotation angles
		double deltaPhi_lab = 2*TMath::Pi()*(this_anal.gen->Rndm());
		double deltaPhi_q = 2*TMath::Pi()*(this_anal.gen->Rndm()); 
		
		//Set input momenta
		TVector3 e_mom = this_event.e4.Vect();
		TVector3 pi_mom[2] = {this_event.pi_1.getPi_q().Vect(), this_event.pi_2.getPi_q().Vect()};
		bool isGoodPion[2] = {this_event.isGoodPion_1, this_event.isGoodPion_2};
		//rotate electron about beam axis and check acceptance
		e_mom.RotateZ(deltaPhi_lab);
		
		bool e_acc = this_anal.anal->checkAcceptance( e_mom.Mag(),rad_to_deg*e_mom.Phi(), rad_to_deg*e_mom.Theta(), 0 ) ;
		if( e_acc < 0 ){continue;}

		for( int i = 0; i < 2; i++ ){
			//rotate pion about q and z
			TVector3 pi_q_mom = pi_mom[i];
			pi_q_mom.RotateZ( deltaPhi_q );
		
			pi_mom[i] = rotate_to_beam_frame( this_event.q, this_event.e4, pi_q_mom );
			pi_mom[i].RotateZ( deltaPhi_lab );
			
			//Check Pion acceptance
			int piType = (int)(i > 0) + 1;
		
			int new_sec = this_anal.anal->checkAcceptance( pi_mom[i].Mag(), rad_to_deg*pi_mom[i].Phi(), rad_to_deg*pi_mom[i].Theta(), piType ) ;
		
			if( new_sec > -1 ) {
				if( this_anal.acc_match ) detectedPion[i] = this_anal.anal->acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
				else detectedPion[i] = true;
			}
			if( new_sec > -1 &&  isGoodPion[i] ){ 				
				accPion[i] = true;
				goodPion[i] =  this_anal.anal->acceptance_match_2d( pi_mom[i].Theta()*rad_to_deg, pi_mom[i].Mag(), new_sec + 1);
			}
			
		}

		if( detectedPion[0] == true && goodPion[0] == true && detectedPion[1] == false ) {
			onePiEvents[0]++;
		}				
		if ( detectedPion[1] == true && goodPion[1] == true && detectedPion[0] == false ) {
			onePiEvents[1]++;
		}
		if( ( detectedPion[0] == true && detectedPion[1] == true ) &&
			( goodPion[0] == true || goodPion[1] == true  ) ){
				twoPiEvents++;
		}

		//check uncertainty
		if( onePiEvents[0] != 0 && twoPiEvents != 0 ){
			corr_err[0] =  sqrt( 1./onePiEvents[0] + 1./twoPiEvents );
		}
		if( onePiEvents[1] != 0 && twoPiEvents != 0 ){
			corr_err[1] =  sqrt( 1./onePiEvents[1] + 1./twoPiEvents );
		}
		//cout<<"CURRENT UNCERTAINTY : "<<corr_err<<endl;	

	}	

	std::array<double, 6> output = {onePiEvents[0], onePiEvents[1], twoPiEvents, corr_err[0], corr_err[1], trials};

	return output;

}



TVector3 rotate_to_beam_frame( TLorentzVector q, TLorentzVector p_e, TVector3 pi_q ){
	TVector3 v = pi_q;
    
	p_e.RotateZ(-q.Phi());
	p_e.RotateY(-q.Theta());

	v.RotateZ( p_e.Phi() );
    v.RotateY( q.Theta() );

    v.RotateZ( q.Phi()  );

	return v;
}

void fill_chunk( TTree * inTree, TTree** tree_chunk, int chunkSize, int i_thread){
	int curr_chunk = -1;
	for( int event = i_thread*chunkSize; event < (i_thread+1)*chunkSize; event++ ){
		if( event >= inTree->GetEntries())break;
		if( event%(chunkSize) == 0 || curr_chunk < 0){
			curr_chunk++;
			std::cout<<"Filling chunk #"<<curr_chunk<<std::endl;
		}
		if( event%(10000) == 0  ) cout<<"Processing Event : "<<event<<" / "<<inTree->GetEntries()<<std::endl;
		inTree->GetEntry(event);
		tree_chunk[i_thread]->Fill();
	}
	tree_chunk[i_thread]->SetDirectory(0);
}


void rotate_chunk( analData this_anal, int i_thread){
	//Set output tree
	auto start = std::chrono::high_resolution_clock::now();

	TFile * inFile = new TFile(this_anal.in_name);
	TFile * f_out = new TFile(Form(this_anal.out_name + "_chunk_%i.root", i_thread), "RECREATE");
	TTree * inChunk = (TTree*)inFile->Get("ePi");
	TTree * outChunk = new TTree("ePi", "(e,e'pi) event  information");

	TLorentzVector beam_out;	
	electron e_out;
	std::vector<pion> pi_out;
	std::vector<bool> goodPiOut;

	double Mx_2pi_out;
	double M_rho_out;
	double rhoWeight[2];
	double corr_err[2];
	int trials;

	electron * e = nullptr;
	TLorentzVector * beam = nullptr;
	std::vector<pion>* pi = nullptr;
	std::vector<bool>* isGoodPion = nullptr;
	std::vector<bool>* isGoodPion_acc = nullptr;
	double Mx_2pi;
	double M_rho;

	inChunk->SetBranchAddress("e", &e);
	inChunk->SetBranchAddress("Mx_2pi", &Mx_2pi);
	inChunk->SetBranchAddress("M_rho", &M_rho);
	inChunk->SetBranchAddress("pi", &pi);
	inChunk->SetBranchAddress("isGoodPion_no_acc", &isGoodPion);
	inChunk->SetBranchAddress("isGoodPion", &isGoodPion_acc);
	inChunk->SetBranchAddress("beam", &beam);


	outChunk->Branch("beam", &beam_out);
	outChunk->Branch("e", &e_out);
	outChunk->Branch("pi", &pi_out);
	outChunk->Branch("isGoodPion", &goodPiOut);
	outChunk->Branch("Mx_2pi", &Mx_2pi_out);
	outChunk->Branch("M_rho", &M_rho_out);

	outChunk->Branch("rhoWeight", rhoWeight, "rhoWeight[2]/D");
	outChunk->Branch("rhoErr", corr_err, "rhoErr[2]/D");

	outChunk->Branch("trials", &trials);

	//std::ofstream txtFile;
	//txtFile.open( "../plotting/rho_data_points.txt" );



	for( int event_count = i_thread*this_anal.chunkSize; event_count < (i_thread+1)*this_anal.chunkSize; event_count++) {
		
		//if(((double)event_count/(double)inChunk->GetEntries())*100 == (int)( (double)event_count/(double)inChunk->GetEntries())*100){
		//	cout<<"Thread "<<i_thread<<" : "<<((double)event_count/(double)inChunk->GetEntries())*100<<std::endl;
		//}

		if((event_count%100) == (0)){
			cout<<"Thread "<<i_thread<<" : "<<event_count<<" / "<<(i_thread+1)*this_anal.chunkSize;
			auto finish = std::chrono::high_resolution_clock::now();
    		std::chrono::duration<double> elapsed = finish - start;
			cout<<" in "<<elapsed.count()/60.<<" minutes\n";
		}		
		if( event_count >= inChunk->GetEntries()) break;
		inChunk->GetEntry(event_count);

		//initialize event variables

		double onePiEvents[2] = {0};
		double twoPiEvents = 0;
		trials = 0;

		pi_out.clear();
		goodPiOut.clear();
		e_out.Clear();

		beam_out = (TLorentzVector) (*beam);
		e_out = (electron)(*e);
		pi_out.push_back((pion)(pi->at(0)));
		pi_out.push_back((pion)(pi->at(1)));
		
		goodPiOut.push_back((bool)isGoodPion->at(0));
		goodPiOut.push_back((bool)isGoodPion->at(1));
		
		Mx_2pi_out = (double)(Mx_2pi);
		M_rho_out = (double)(M_rho);
		rhoWeight[0] = 0; // (pi+, pi-)
		rhoWeight[1] = 0; // (pi+, pi-)                     
		corr_err[0] = 999;
		corr_err[1] = 999;


		//Look at particles in acceptance map
		if( this_anal.anal->checkAcceptance( e->get3Momentum().Mag(),
					rad_to_deg*e->get3Momentum().Phi(), 
					rad_to_deg*e->get3Momentum().Theta(), 0 )  < 0){continue;}
		if( this_anal.anal->checkAcceptance( (*pi)[0].get3Momentum().Mag(), 
					rad_to_deg*(*pi)[0].get3Momentum().Phi(), 
					rad_to_deg*(*pi)[0].get3Momentum().Theta(), 
					(int)( (*pi)[0].getCharge() < 0 ) + 1 ) < 0 ) { continue; }
		if( this_anal.anal->checkAcceptance( (*pi)[1].get3Momentum().Mag(), 
					rad_to_deg*(*pi)[1].get3Momentum().Phi(), 
					rad_to_deg*(*pi)[1].get3Momentum().Theta(), 
					(int)( (*pi)[1].getCharge() < 0 ) + 1 ) < 0  ) { continue; }
		//Restrict ROI			
		if( this_anal.acc_match && !this_anal.anal->applyAcceptanceMatching((*pi)[0], 2) ){continue;}
		if( this_anal.acc_match && !this_anal.anal->applyAcceptanceMatching((*pi)[1], 2) ){continue;}
		//if( !isGoodPion_acc[0] && !isGoodPion_acc[1] ){continue;} //If event wouldn't be in our final sample, continue 
		if( Mx_2pi < 0 || Mx_2pi > 2.5 ){continue;}
		if( M_rho < 0 ){continue;}
		
		eventData this_event = {e->get4Momentum(), e->getQ(), (*pi)[0], (*pi)[1], (*isGoodPion)[0], (*isGoodPion)[1]};
		std::array<double, 6> rotOut = rotateEvent(this_event, this_anal );
		
		
		onePiEvents[0] = rotOut[0];
		onePiEvents[1] = rotOut[1];
		twoPiEvents = rotOut[2];
		corr_err[0] = rotOut[3];
		corr_err[1] = rotOut[4];
		trials = rotOut[5];

		//Check for weird guys
		for( int i = 0; i < 2; i++ ){
			if( onePiEvents[i] == 0 || twoPiEvents == 0 ){
				if( (*isGoodPion_acc)[i] ){
					rhoWeight[i] = 1.;
				}
				else{ rhoWeight[i] = 0.;}
			}
			else if( (*isGoodPion_acc)[i] ){
				rhoWeight[i] =  1. + onePiEvents[i]/twoPiEvents ;
			}
			else{
				rhoWeight[i] = onePiEvents[i]/twoPiEvents ;
			}
		}
		outChunk->Fill();

	}
	
	f_out->cd();
	outChunk->Write();
	f_out->Close();


	
}

