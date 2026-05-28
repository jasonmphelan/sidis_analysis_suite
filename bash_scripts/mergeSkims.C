// mergeSkims.C
// Merges two GEMC final_skim trees where one (gen_skimmer output) is missing
// branches that the other (test_skimmer output) has.
// Missing branches are added and filled with safe defaults before hadd.
//
// Usage:
//   root -l -b -q 'mergeSkims.C("output.root", "test_skim.root", "gen_skim.root")'
//   root -l -b -q 'mergeSkims.C("output.root", "test_skim.root", "gen_skim.root", 1)'
//   (last arg = default value for "target" branch, default 0)

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <vector>

using std::cout;
using std::cerr;

void patchMissingBranches( TString fileName, int defaultTarget ){
    cout << "Patching missing branches in " << fileName << "\n";
    TFile* f = TFile::Open( fileName, "UPDATE" );
    if( !f || f->IsZombie() ){ cerr << "Cannot open " << fileName << "\n"; return; }
    TTree* t = (TTree*) f->Get("ePi");
    if( !t ){ cerr << "No ePi tree in " << fileName << "\n"; return; }

    Long64_t n = t->GetEntries();

    // Scalar branches
    int    target    = defaultTarget;
    int    accCharge = 0;
    bool   isDIS     = false;
    int    Nka = 0, Nkaps = 0, Nakms = 0;
    int    Npr = 0, Nprps = 0, Naprms = 0;

    // Vector branches (empty per event)
    std::vector<int>*    pidRICH  = new std::vector<int>();
    std::vector<double>* probRICH = new std::vector<double>();
    std::vector<double>* nRICH    = new std::vector<double>();

    // Add only branches that are actually missing
    auto addInt = [&]( const char* name, int* addr ){
        if( !t->GetBranch(name) ){
            cout << "  Adding branch: " << name << "\n";
            TBranch* b = t->Branch(name, addr, Form("%s/I", name));
            for( Long64_t i = 0; i < n; i++ ) b->Fill();
        }
    };
    auto addVecI = [&]( const char* name, std::vector<int>** addr ){
        if( !t->GetBranch(name) ){
            cout << "  Adding branch: " << name << "\n";
            TBranch* b = t->Branch(name, addr);
            for( Long64_t i = 0; i < n; i++ ) b->Fill();
        }
    };
    auto addVecD = [&]( const char* name, std::vector<double>** addr ){
        if( !t->GetBranch(name) ){
            cout << "  Adding branch: " << name << "\n";
            TBranch* b = t->Branch(name, addr);
            for( Long64_t i = 0; i < n; i++ ) b->Fill();
        }
    };

    addInt("target",   &target);
    addInt("accCharge",&accCharge);
    // isDIS must be Bool_t to match test_skimmer output.
    // If missing: add it. If present as Int_t (from a previous bad patch): rewrite tree.
    TBranch* isDISbranch = t->GetBranch("isDIS");
    bool isDIS_wrongType = isDISbranch &&
        TString(isDISbranch->GetLeaf("isDIS")->GetTypeName()) != "Bool_t";

    if( isDIS_wrongType ){
        cout << "  Fixing isDIS type (Int_t -> Bool_t) — requires full tree rewrite\n";

        // Read old Int_t isDIS values
        int isDIS_old = 0;
        t->SetBranchAddress("isDIS", &isDIS_old);

        // Clone structure without isDIS into a temp tree
        t->SetBranchStatus("isDIS", 0);
        TTree* tNew = t->CloneTree(0);
        t->SetBranchStatus("isDIS", 1);

        // Add isDIS back as Bool_t
        TBranch* bNew = tNew->Branch("isDIS", &isDIS, "isDIS/O");

        for( Long64_t i = 0; i < n; i++ ){
            t->GetEntry(i);
            isDIS = (bool) isDIS_old;
            tNew->Fill();
        }

        f->cd();
        tNew->Write("ePi_fixed", TObject::kOverwrite);
        f->Close();

        // Replace the original tree with the fixed one via rename
        TFile* f2 = TFile::Open( fileName, "UPDATE" );
        f2->Delete("ePi;*");
        TTree* tFixed = (TTree*) f2->Get("ePi_fixed");
        tFixed->SetName("ePi");
        tFixed->Write("", TObject::kOverwrite);
        f2->Delete("ePi_fixed;*");
        f2->Close();
        cout << "  Done fixing isDIS type\n";
        return;  // already closed, skip bottom Close()
    } else if( !isDISbranch ){
        cout << "  Adding branch: isDIS\n";
        TBranch* b = t->Branch("isDIS", &isDIS, "isDIS/O");
        for( Long64_t i = 0; i < n; i++ ) b->Fill();
    }
    addInt("Nka",      &Nka);
    addInt("Nkaps",    &Nkaps);
    addInt("Nakms",    &Nakms);
    addInt("Npr",      &Npr);
    addInt("Nprps",    &Nprps);
    addInt("Naprms",   &Naprms);
    addVecI("pidRICH",  &pidRICH);
    addVecD("probRICH", &probRICH);
    addVecD("nRICH",    &nRICH);

    t->Write("", TObject::kOverwrite);
    f->Close();
    cout << "Done patching " << fileName << "\n";
}

void mergeSkims( TString outName, TString file1, TString file2, int defaultTarget = 0 ){

    // Inspect both files for missing branches
    TString branchNames[] = {
        "target","accCharge","isDIS",
        "Nka","Nkaps","Nakms","Npr","Nprps","Naprms",
        "pidRICH","probRICH","nRICH"
    };
    int nBranches = sizeof(branchNames)/sizeof(branchNames[0]);

    auto needsPatch = [&]( TString fname ) -> bool {
        TFile* f = TFile::Open(fname, "READ");
        TTree* t = (TTree*) f->Get("ePi");
        bool needs = false;
        // Missing branch?
        for( int i = 0; i < nBranches; i++ )
            if( !t->GetBranch(branchNames[i]) ){ needs = true; break; }
        // isDIS present but wrong type?
        TBranch* b = t->GetBranch("isDIS");
        if( b && TString(b->GetLeaf("isDIS")->GetTypeName()) != "Bool_t" )
            needs = true;
        f->Close();
        return needs;
    };

    if( needsPatch(file1) ) patchMissingBranches( file1, defaultTarget );
    if( needsPatch(file2) ) patchMissingBranches( file2, defaultTarget );

    TString cmd = "hadd -f " + outName + " " + file1 + " " + file2;
    cout << "Running: " << cmd << "\n";
    gSystem->Exec( cmd );
    cout << "Output: " << outName << "\n";
}
