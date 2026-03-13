// Reads output of final_skimmer and rewrites isGoodPion branches with tighter
// cuts: full kinematic cuts + detector cuts (fiducials, chi2, vertex) + acceptance.

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TEventList.h"
#include "cut_values.h"
#include "electron.h"
#include "pion.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"

using std::cerr;
using std::cout;
using namespace cutVals;

int main(int argc, char** argv) {

	auto start = std::chrono::high_resolution_clock::now();

	if (argc < 6) {
		cerr << "Incorrect number of arguments. Please use:\n";
		cerr << "./tight_reskim [Input File] [Output File] [Beam Energy] [Target] [RunType]\n";
		cerr << "       Target:  0 = RGB/deuterium, 1 = RGA/proton\n";
		cerr << "       RunType: 0 = data, 1 = MC\n";
		return -1;
	}

	TString in_name  = argv[1];
	TString out_name = argv[2];
	double  EBeam    = atof(argv[3]);
	int     target   = atoi(argv[4]);
	int     runType  = atoi(argv[5]);

	cerr << "Input: " << in_name << "  Output: " << out_name
	     << "  EBeam: " << EBeam << "  Target: " << target
	     << "  RunType: " << runType << "\n";

	TFile* inFile = new TFile(in_name, "READ");
	if (!inFile || inFile->IsZombie()) {
		cerr << "Failed to open input file: " << in_name << "\n";
		return -1;
	}

	TTree* inTree = (TTree*) inFile->Get("ePi");
	if (!inTree) {
		cerr << "Could not find tree 'ePi' in " << in_name << "\n";
		return -1;
	}

	analyzer anal(0, -1);
	anal.setAnalyzerLevel(0);
	anal.setTarget(target);
	anal.loadMatchingFunctions(target == 1 ? "matchCut2D_rga.root" : "matchCut2D_map.root");
	anal.loadMatchingFunctions3D();
	anal.loadAcceptanceMapContinuous((TString)_DATA + (TString)"/acceptance_map/acceptanceMap_allE_final.root");

	// Set up input branch addresses
	electron*             e      = nullptr;
	std::vector<pion>*    pi     = nullptr;
	genElectron*          e_gen  = nullptr;
	std::vector<genPion>* pi_gen = nullptr;

	inTree->SetBranchAddress("e",  &e);
	inTree->SetBranchAddress("pi", &pi);
	if (runType == 1) {
		inTree->SetBranchAddress("e_gen",  &e_gen);
		inTree->SetBranchAddress("pi_gen", &pi_gen);
	}

	// Clone tree structure but exclude the old isGoodPion branches; we will
	// recompute them and add fresh copies to the output.
	inTree->SetBranchStatus("isGoodPion",        0);
	inTree->SetBranchStatus("isGoodPion_no_acc", 0);
	inTree->SetBranchStatus("isGoodPion_3d",     0);
	inTree->SetBranchStatus("isGoodGenPion",     0);

	TFile* outFile = new TFile(out_name, "RECREATE");
	TTree* outTree = inTree->CloneTree(0);

	// Re-enable so GetEntry reads them if needed in future, though we don't use them here
	inTree->SetBranchStatus("isGoodPion",        1);
	inTree->SetBranchStatus("isGoodPion_no_acc", 1);
	inTree->SetBranchStatus("isGoodPion_3d",     1);
	inTree->SetBranchStatus("isGoodGenPion",     1);

	std::vector<bool> isGoodPion;
	std::vector<bool> isGoodPion_no_acc;
	std::vector<bool> isGoodPion_3d;
	std::vector<bool> isGoodGenPion;

	outTree->Branch("isGoodPion",        &isGoodPion);
	outTree->Branch("isGoodGenPion",     &isGoodGenPion);
	outTree->Branch("isGoodPion_no_acc", &isGoodPion_no_acc);
	outTree->Branch("isGoodPion_3d",     &isGoodPion_3d);

	long long nEvents = inTree->GetEntries();
	cout << "Starting event loop over " << nEvents << " events\n";

	for (long long ev = 0; ev < nEvents; ev++) {

		if (ev % 100000 == 0) cout << ev << " / " << nEvents << "\n";

		isGoodPion.clear();
		isGoodPion_no_acc.clear();
		isGoodPion_3d.clear();
		isGoodGenPion.clear();

		inTree->GetEntry(ev);

		// Tighter electron: kinematic + full detector cuts
		bool e_pass     = anal.applyElectronKinematicCuts(*e) &&
		                  anal.applyElectronDetectorCuts(*e);
		bool e_gen_pass = false;
		if (runType == 1) e_gen_pass = anal.applyElectronKinematicCuts(*e_gen);

		if (!e_pass && !e_gen_pass) continue;

		bool hasGoodReco = false;
		bool hasGoodGen  = false;

		for (int i = 0; i < (int) pi->size(); i++) {

			isGoodPion.push_back(false);
			isGoodPion_no_acc.push_back(false);
			isGoodPion_3d.push_back(false);
			isGoodGenPion.push_back(false);

			int chargeIdx = (int)((*pi)[i].getCharge() < 1);

			// Tighter pion: kinematic + full detector cuts
			bool pi_pass     = anal.applyPionKinematicCuts((*pi)[i]) &&
			                   anal.applyPionDetectorCuts((*pi)[i], *e);
			bool pi_gen_pass = false;
			if (runType == 1) pi_gen_pass = anal.applyPionKinematicCuts((*pi_gen)[i]);

			if (!pi_pass && !pi_gen_pass) continue;

			isGoodGenPion[i]     = (pi_gen_pass && e_gen_pass);
			isGoodPion_no_acc[i] = (pi_pass && e_pass);

			if (isGoodGenPion[i])     hasGoodGen  = true;
			if (isGoodPion_no_acc[i]) hasGoodReco = true;

			double p_pi     = (*pi)[i].get3Momentum().Mag();
			double theta_pi = (*pi)[i].get3Momentum().Theta() * rad_to_deg;
			double phi_pi   = (*pi)[i].get3Momentum().Phi()   * rad_to_deg;

			if (anal.applyAcceptanceMap(e->get3Momentum().Mag(),
			                            rad_to_deg * e->get3Momentum().Phi(),
			                            rad_to_deg * e->get3Momentum().Theta(), 0) < 0) continue;
			if (anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, chargeIdx + 1) < 0) continue;

			if (anal.applyAcceptanceMatching((*pi)[i], 2) && (pi_pass && e_pass))
				isGoodPion[i] = true;

			if (anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, 1) >= 0 &&
			    anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, 2) >= 0 &&
			    (pi_pass && e_pass))
				isGoodPion_3d[i] = true;
		}

		if (hasGoodReco || hasGoodGen)
			outTree->Fill();
	}

	cout << "Writing to file\n";
	outFile->cd();
	outTree->Write();
	outFile->Close();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Done. Elapsed Time: " << elapsed.count() << "\n";

	return 0;
}
