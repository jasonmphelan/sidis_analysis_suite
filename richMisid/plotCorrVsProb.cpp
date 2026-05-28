/**
 * plotCorrVsProb.cpp
 *
 * Reads the ePi TTree from a kaon skim file and plots the RICH kaon
 * misidentification correction C(p,θ) — loaded from the Stage-2 output —
 * as a function of probICH and momentum.
 *
 * For each good pion the correction factor is looked up from
 * hCorrFactor_{Kp,Km}(p, θ) and filled into:
 *   • TH2F  hCorr_vs_prob_p  – C vs (probICH, p)         [2D]
 *   • TProfile hCorr_prof_prob – mean C vs probICH        [profile]
 *   • TProfile hCorr_prof_p    – mean C vs p              [profile]
 *   • TH1F  hProb              – probICH distribution
 *   • TH2F  hPID_prob          – pidRICH vs probICH
 *
 * Usage:
 *   ./plotCorrVsProb <tree.root> <corrections.root> [output.root]
 *
 * The corrections file is the Stage-2 output from richKaonMisidCorr.
 * Histograms are written to output.root (default: corrVsProb.root) and
 * canvases are saved as PDFs alongside it.
 */

#include <iostream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "pion.h"
#include "electron.h"
#include "constants.h"

using namespace constants;

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "\nUsage: " << argv[0]
                  << " <tree.root> <corrections.root> [output.root]\n\n";
        return -1;
    }

    const char* tree_name = argv[1];
    const char* corr_name = argv[2];
    const char* out_name  = (argc > 3) ? argv[3] : "corrVsProb.root";

    std::cerr << "Tree file   : " << tree_name << "\n";
    std::cerr << "Corrections : " << corr_name << "\n";
    std::cerr << "Output      : " << out_name  << "\n";

    // ── Open correction file and load hCorrFactor histograms ─────────────────
    TFile* corrFile = TFile::Open(corr_name);
    if (!corrFile || corrFile->IsZombie()) {
        std::cerr << "ERROR: cannot open corrections file: " << corr_name << "\n";
        return -1;
    }

    static const char* cstr[2] = {"Kp", "Km"};
    TH2F* hCorr[2];
    for (int ic = 0; ic < 2; ++ic) {
        hCorr[ic] = dynamic_cast<TH2F*>(
            corrFile->Get(Form("hCorrFactor_%s", cstr[ic])));
        if (!hCorr[ic]) {
            std::cerr << "ERROR: hCorrFactor_" << cstr[ic]
                      << " not found in " << corr_name << "\n";
            return -1;
        }
        hCorr[ic]->SetDirectory(nullptr);
    }
    corrFile->Close();
    std::cerr << "Correction histograms loaded.\n";

    // ── Open tree ─────────────────────────────────────────────────────────────
    TFile* treeFile = TFile::Open(tree_name);
    if (!treeFile || treeFile->IsZombie()) {
        std::cerr << "ERROR: cannot open tree file: " << tree_name << "\n";
        return -1;
    }
    TTree* tree = dynamic_cast<TTree*>(treeFile->Get("ePi"));
    if (!tree) {
        std::cerr << "ERROR: TTree 'ePi' not found in " << tree_name << "\n";
        return -1;
    }
    std::cerr << "Tree entries: " << tree->GetEntries() << "\n";

    // ── Set up branches ───────────────────────────────────────────────────────
    TTreeReader reader(tree);

    TTreeReaderValue<electron>           e      (reader, "e");
    TTreeReaderArray<pion>               pi     (reader, "pi");
    TTreeReaderArray<bool>               isGood (reader, "isGoodPion_no_acc");
    TTreeReaderArray<int>                pidRICH(reader, "pidRICH");
    TTreeReaderArray<double>             probICH(reader, "probICH");

    // ── Declare output histograms (per charge) ────────────────────────────────
    // C vs (probICH, p)
    TH2F* hCorr_vs_prob_p[2];
    // Mean C as profiles
    TProfile* hCorr_prof_prob[2];
    TProfile* hCorr_prof_p[2];
    // probICH distributions
    TH1F* hProb[2];
    // pidRICH vs probICH
    TH2F* hPID_prob[2];

    static const char* ctex[2] = {"K^{+}", "K^{-}"};

    for (int ic = 0; ic < 2; ++ic) {
        hCorr_vs_prob_p[ic] = new TH2F(
            Form("hCorr_vs_prob_p_%s", cstr[ic]),
            Form("C(p,#theta) vs prob_{ICH} and p  (%s);"
                 "prob_{ICH};p (GeV/c);C(p,#theta)", ctex[ic]),
            50, 0.0, 1.0,   // probICH axis
            30, 1.25, 5.0); // p axis
        hCorr_vs_prob_p[ic]->Sumw2();

        hCorr_prof_prob[ic] = new TProfile(
            Form("hCorr_prof_prob_%s", cstr[ic]),
            Form("Mean C(p,#theta) vs prob_{ICH}  (%s);"
                 "prob_{ICH};<C(p,#theta)>", ctex[ic]),
            50, 0.0, 1.0);

        hCorr_prof_p[ic] = new TProfile(
            Form("hCorr_prof_p_%s", cstr[ic]),
            Form("Mean C(p,#theta) vs p  (%s);"
                 "p (GeV/c);<C(p,#theta)>", ctex[ic]),
            30, 1.25, 5.0);

        hProb[ic] = new TH1F(
            Form("hProb_%s", cstr[ic]),
            Form("prob_{ICH} distribution  (%s);"
                 "prob_{ICH};Counts", ctex[ic]),
            50, 0.0, 1.0);

        hPID_prob[ic] = new TH2F(
            Form("hPID_prob_%s", cstr[ic]),
            Form("pid_{RICH} vs prob_{ICH}  (%s);"
                 "prob_{ICH};pid_{RICH}", ctex[ic]),
            50, 0.0, 1.0,
            30, -15.5, 14.5);
    }

    // ── Event loop ────────────────────────────────────────────────────────────
    long long nProc = 0, nFilled = 0;

    while (reader.Next()) {
        if (nProc % 200000 == 0)
            std::cout << "Events processed: " << nProc << "\n";
        ++nProc;

        int nPions = static_cast<int>(pi.end() - pi.begin());
        for (int i = 0; i < nPions; ++i) {

            if (!isGood[i]) continue;

            double p_pi    = pi[i].get3Momentum().Mag();
            double theta_pi = pi[i].get3Momentum().Theta() * rad_to_deg;
            int    ic      = (pi[i].getCharge() < 0) ? 1 : 0;

            // Look up correction from hCorrFactor(p, θ)
            int binX = hCorr[ic]->GetXaxis()->FindBin(p_pi);
            int binY = hCorr[ic]->GetYaxis()->FindBin(theta_pi);
            // Clamp to valid range
            binX = std::max(1, std::min(binX, hCorr[ic]->GetNbinsX()));
            binY = std::max(1, std::min(binY, hCorr[ic]->GetNbinsY()));
            double corr = hCorr[ic]->GetBinContent(binX, binY);
            if (corr <= 0.0) continue;  // bin not filled / outside acceptance

            double prob = probICH[i];
            int    pid  = pidRICH[i];

            hCorr_vs_prob_p[ic]->Fill(prob, p_pi, corr);
            hCorr_prof_prob [ic]->Fill(prob, corr);
            hCorr_prof_p    [ic]->Fill(p_pi, corr);
            hProb           [ic]->Fill(prob);
            hPID_prob       [ic]->Fill(prob, pid);

            ++nFilled;
        }
    }

    std::cout << "\nDone.  Events: " << nProc
              << "  Pions filled: " << nFilled << "\n";

    treeFile->Close();

    // ── Draw and save ─────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    TFile* outFile = new TFile(out_name, "RECREATE");

    for (int ic = 0; ic < 2; ++ic) {

        // 1) C vs (probICH, p) — 2D
        {
            TCanvas* c1 = new TCanvas(
                Form("cCorr_vs_prob_p_%s", cstr[ic]), "", 1600, 700);
            c1->SetLeftMargin(0.10); c1->SetRightMargin(0.14);
            c1->SetBottomMargin(0.13); c1->SetTopMargin(0.10);
            hCorr_vs_prob_p[ic]->Draw("COLZ");
            c1->SaveAs(Form("corrVsProb_2d_%s.pdf", cstr[ic]));
            outFile->cd(); c1->Write(); delete c1;
        }

        // 2) Profile: mean C vs probICH
        {
            TCanvas* c2 = new TCanvas(
                Form("cCorr_prof_prob_%s", cstr[ic]), "", 900, 650);
            c2->SetLeftMargin(0.13); c2->SetBottomMargin(0.13);
            hCorr_prof_prob[ic]->SetLineColor(kBlue + 1);
            hCorr_prof_prob[ic]->SetLineWidth(2);
            hCorr_prof_prob[ic]->SetMarkerStyle(20);
            hCorr_prof_prob[ic]->SetMarkerColor(kBlue + 1);
            hCorr_prof_prob[ic]->Draw("E");
            c2->SaveAs(Form("corrVsProb_profProb_%s.pdf", cstr[ic]));
            outFile->cd(); c2->Write(); delete c2;
        }

        // 3) Profile: mean C vs p
        {
            TCanvas* c3 = new TCanvas(
                Form("cCorr_prof_p_%s", cstr[ic]), "", 900, 650);
            c3->SetLeftMargin(0.13); c3->SetBottomMargin(0.13);
            hCorr_prof_p[ic]->SetLineColor(kRed + 1);
            hCorr_prof_p[ic]->SetLineWidth(2);
            hCorr_prof_p[ic]->SetMarkerStyle(20);
            hCorr_prof_p[ic]->SetMarkerColor(kRed + 1);
            hCorr_prof_p[ic]->Draw("E");
            c3->SaveAs(Form("corrVsProb_profP_%s.pdf", cstr[ic]));
            outFile->cd(); c3->Write(); delete c3;
        }

        // 4) probICH distribution
        {
            TCanvas* c4 = new TCanvas(
                Form("cProb_%s", cstr[ic]), "", 900, 650);
            c4->SetLeftMargin(0.13); c4->SetBottomMargin(0.13);
            hProb[ic]->SetLineColor(kBlue + 1);
            hProb[ic]->SetLineWidth(2);
            hProb[ic]->Draw("HIST");
            c4->SaveAs(Form("corrVsProb_probDist_%s.pdf", cstr[ic]));
            outFile->cd(); c4->Write(); delete c4;
        }

        // 5) pidRICH vs probICH
        {
            TCanvas* c5 = new TCanvas(
                Form("cPID_prob_%s", cstr[ic]), "", 1200, 700);
            c5->SetLeftMargin(0.10); c5->SetRightMargin(0.14);
            c5->SetBottomMargin(0.13); c5->SetTopMargin(0.10);
            hPID_prob[ic]->Draw("COLZ");
            c5->SaveAs(Form("corrVsProb_pidVsProb_%s.pdf", cstr[ic]));
            outFile->cd(); c5->Write(); delete c5;
        }

        outFile->cd();
        hCorr_vs_prob_p [ic]->Write();
        hCorr_prof_prob [ic]->Write();
        hCorr_prof_p    [ic]->Write();
        hProb           [ic]->Write();
        hPID_prob       [ic]->Write();
    }

    outFile->Close();
    std::cerr << "Output written to: " << out_name << "\n";
    return 0;
}
