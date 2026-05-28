/**
 * pidRICHHists.cpp
 *
 * Loops over the ePi TTree and fills, for each probICH cut value:
 *
 *   hAll_{cut}   – (p, θ) of all good K+ with probICH > cut
 *   hPi_{cut}    – (p, θ) of tracks where pidRICH == -211 (pion-tagged)
 *   hRatio_{cut} – hPi / hAll  (pion contamination fraction)
 *
 * One directory per cut value is written to the output ROOT file.
 * A multi-page PDF is produced with one page per cut.
 *
 * Usage:
 *   ./pidRICHHists <tree.root> [output.root]
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "pion.h"
#include "electron.h"
#include "constants.h"
#include "cut_values.h"
#include "analyzer.h"

using namespace constants;

// ── Histogram binning ─────────────────────────────────────────────────────────
static constexpr int    N_P_BINS  = 50;
static constexpr double P_MIN     = 1.25;
static constexpr double P_MAX     = 5.00;

static constexpr int    N_TH_BINS = 75;
static constexpr double TH_MIN    =  5.0;
static constexpr double TH_MAX    = 35.0;

// ── probICH cut values to scan ────────────────────────────────────────────────
static const std::vector<double> PROB_CUTS = { -.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "\nUsage: " << argv[0]
                  << " <tree.root> [output.root]\n\n";
        return -1;
    }

    const char* tree_name = argv[1];
    const char* out_name  = (argc > 2) ? argv[2] : "pidRICHHists.root";

    std::cerr << "Tree file : " << tree_name << "\n";
    std::cerr << "Output    : " << out_name  << "\n";

    // ── Acceptance map ────────────────────────────────────────────────────────
    analyzer anal(0, -1);
    anal.setAnalyzerLevel(0);
    anal.loadAcceptanceMapContinuous(
        (TString)_DATA + "/acceptance_map/acceptanceMap_allE_final.root");

    // ── Open tree ─────────────────────────────────────────────────────────────
    TChain* chain = new TChain("ePi");
    chain->Add(tree_name);
    long long nEntries = chain->GetEntries();
    std::cerr << "Entries   : " << nEntries << "\n";
    if (nEntries == 0) { std::cerr << "ERROR: no entries.\n"; return -1; }

    TTreeReader reader(chain);
    TTreeReaderValue<electron>       e      (reader, "e");
    TTreeReaderArray<pion>           pi     (reader, "pi");
    TTreeReaderArray<bool>           isGood (reader, "isGoodPion_no_acc");
    TTreeReaderArray<int>            pidRICH(reader, "pidRICH");
    TTreeReaderArray<double>         probICH(reader, "probRICH");

    // ── Allocate histograms for each cut ──────────────────────────────────────
    int nCuts = static_cast<int>(PROB_CUTS.size());
    std::vector<TH2F*> hAll(nCuts), hPi(nCuts), hRatio(nCuts);

    for (int ic = 0; ic < nCuts; ++ic) {
        double cut = PROB_CUTS[ic];
        std::string tag = Form("prob%.2g", cut);
        for (char& ch : tag) if (ch == '.') ch = 'p';

        hAll[ic] = new TH2F(
            Form("hAll_%s",   tag.c_str()),
            Form("All K^{+}  probICH > %.2g;p (GeV/c);#theta (deg)", cut),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);

        hPi[ic] = new TH2F(
            Form("hPi_%s",    tag.c_str()),
            Form("pid==-211  probICH > %.2g;p (GeV/c);#theta (deg)", cut),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);

        hRatio[ic] = new TH2F(
            Form("hRatio_%s", tag.c_str()),
            Form("Pion fraction  probICH > %.2g;p (GeV/c);#theta (deg)", cut),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);
        hRatio[ic]->Sumw2();
    }

    // ── Event loop ────────────────────────────────────────────────────────────
    long long nProc = 0;
    while (reader.Next()) {
        if (nProc % 200000 == 0)
            std::cout << "Events processed: " << nProc << "\n";
        ++nProc;

        // Electron acceptance check
        double p_e     = e->get3Momentum().Mag();
        double phi_e   = e->get3Momentum().Phi()   * rad_to_deg;
        double theta_e = e->get3Momentum().Theta() * rad_to_deg;
        if (anal.applyAcceptanceMap(p_e, phi_e, theta_e, 0) < 0) continue;

        int nPions = static_cast<int>(pi.end() - pi.begin());
        for (int i = 0; i < nPions; ++i) {

            if (!isGood[i]) continue;
            if (pi[i].getCharge() < 0) continue;  // K+ only

            double p_pi     = pi[i].get3Momentum().Mag();
            double theta_pi = pi[i].get3Momentum().Theta() * rad_to_deg;
            double phi_pi   = pi[i].get3Momentum().Phi()   * rad_to_deg;

            if (p_pi    < P_MIN  || p_pi    >= P_MAX)  continue;
            if (theta_pi < TH_MIN || theta_pi >= TH_MAX) continue;

            // Pion acceptance check (chargeIdx 1 = pi+)
            if (anal.applyAcceptanceMap(p_pi, phi_pi, theta_pi, 1) < 0) continue;
            if( pi[i].getBeta_rich() < .01 || pidRICH[i] == 1) continue;
            double prob = probICH[i];
            int    pid  = pidRICH[i];

            for (int ic = 0; ic < nCuts; ++ic) {
                if (prob <= PROB_CUTS[ic]) continue;
                hAll[ic]->Fill(p_pi, theta_pi);
                if (pid == 211)
                    hPi[ic]->Fill(p_pi, theta_pi);
            }
        }
    }

    std::cout << "\nDone.  Events processed: " << nProc << "\n";

    // ── Compute ratios ────────────────────────────────────────────────────────
    for (int ic = 0; ic < nCuts; ++ic) {
        hRatio[ic]->Add(hPi[ic]);
        hRatio[ic]->Divide(hAll[ic]);
    }

    // ── Draw multi-page PDF ───────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TString pdfName = TString(out_name).ReplaceAll(".root", ".pdf");
    TCanvas* c = new TCanvas("c", "", 1800, 650);
    c->Print(pdfName + "[");

    for (int ic = 0; ic < nCuts; ++ic) {
        c->Clear();
        c->Divide(3, 1, 0.003, 0.003);

        auto stylepad = [](TVirtualPad* p) {
            p->SetLeftMargin(0.12); p->SetRightMargin(0.16);
            p->SetBottomMargin(0.13); p->SetTopMargin(0.12);
        };

        c->cd(1); stylepad(gPad); hAll  [ic]->Draw("COLZ");
        c->cd(2); stylepad(gPad); hPi   [ic]->Draw("COLZ");
        c->cd(3); stylepad(gPad); hRatio[ic]->Draw("COLZ");

        c->cd(0);
        TLatex lx;
        lx.SetNDC(); lx.SetTextSize(0.030);
        lx.DrawLatex(0.38, 0.965, Form("probICH > %.2g", PROB_CUTS[ic]));

        c->Print(pdfName);
    }

    c->Print(pdfName + "]");
    std::cout << "PDF written to: " << pdfName << "\n";

    // ── Write to ROOT file ────────────────────────────────────────────────────
    TFile* outFile = new TFile(out_name, "RECREATE");
    for (int ic = 0; ic < nCuts; ++ic) {
        std::string tag = Form("prob%.2g", PROB_CUTS[ic]);
        for (char& ch : tag) if (ch == '.') ch = 'p';

        outFile->mkdir(tag.c_str());
        outFile->cd(tag.c_str());
        hAll  [ic]->Write();
        hPi   [ic]->Write();
        hRatio[ic]->Write();
    }

    outFile->Close();
    std::cout << "ROOT written to: " << out_name << "\n";

    delete chain;
    return 0;
}
