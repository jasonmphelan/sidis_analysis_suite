/**
 * drawTree.C
 *
 * ROOT macro: opens a TTree, loops over probRICH cuts, and for each cut
 * produces hPi/hAll (pion-ID ratio vs p,theta).  Results are saved to
 * a multi-page PDF and to a ROOT file (one directory per cut value).
 *
 * Usage (from ROOT prompt):
 *   .x drawTree.C("../trees/kFin.root")
 *
 * Or from the shell:
 *   root -l -q 'drawTree.C("../trees/kFin.root")'
 */

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <string>

void drawTree(const char* fileName = "../trees/final",
              const char* treeName = "ePi",
              const char* outName  = "drawTree_out.root")
{
    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << fileName << "\n";
        return;
    }

    TTree* t = dynamic_cast<TTree*>(f->Get(treeName));
    if (!t) {
        std::cerr << "ERROR: tree '" << treeName << "' not found.\n";
        return;
    }

    std::cout << "Opened '" << treeName << "' with "
              << t->GetEntries() << " entries.\n";

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // ── probRICH cut values to scan ───────────────────────────────────────────
    std::vector<double> probCuts = { 0.0, 0.1, 0.2, 0.3, 0.5 };
    // ─────────────────────────────────────────────────────────────────────────

    TFile* outFile = new TFile(outName, "RECREATE");
    TCanvas* c = new TCanvas("c", "", 900, 650);
    TString pdfName = TString(outName).ReplaceAll(".root", ".pdf");
    c->Print(pdfName + "[");

    for (double cut : probCuts) {

        std::string tag = Form("prob%.2g", cut);
        // Replace '.' with 'p' for ROOT name safety
        for (char& ch : tag) if (ch == '.') ch = 'p';

        TH2F* hAll = new TH2F(Form("hAll_%s", tag.c_str()),
            Form("All K  (probRICH > %.2g);p (GeV/c);#theta (deg)", cut),
            25, 1.25, 5.0, 60, 5, 35);

        TH2F* hPi = new TH2F(Form("hPi_%s", tag.c_str()),
            Form("pid==211  (probRICH > %.2g);p (GeV/c);#theta (deg)", cut),
            25, 1.25, 5.0, 60, 5, 35);

        TH2F* hRatio = new TH2F(Form("hRatio_%s", tag.c_str()),
            Form("Pion fraction  (probRICH > %.2g);p (GeV/c);#theta (deg)", cut),
            25, 1.25, 5.0, 60, 5, 35);

        // Fill hAll and hPi via Draw
        t->Draw(Form("pi.pi3.Theta()*180./3.14:pi.pi3.Mag()>>hAll_%s", tag.c_str()),
                Form("nRICH>2 && probRICH>%g && pi.Charge>0 && e.vt.Z()>-7 && e.vt.Z()<2", cut),
                "goff");

        t->Draw(Form("pi.pi3.Theta()*180./3.14:pi.pi3.Mag()>>hPi_%s", tag.c_str()),
                Form("nRICH>2 && probRICH>%g && pi.Charge>0 && pidRICH==211 && e.vt.Z()>-7 && e.vt.Z()<2", cut),
                "goff");

        // Compute ratio
        hRatio->Add(hPi);
        hRatio->Divide(hAll);

        // ── Draw page: hAll | hPi | ratio ────────────────────────────────────
        c->Clear();
        c->Divide(3, 1, 0.003, 0.003);

        auto stylepad = [](TVirtualPad* p) {
            p->SetLeftMargin(0.12); p->SetRightMargin(0.15);
            p->SetBottomMargin(0.13); p->SetTopMargin(0.12);
        };

        c->cd(1); stylepad(gPad);
        hAll->Draw("COLZ");

        c->cd(2); stylepad(gPad);
        hPi->Draw("COLZ");

        c->cd(3); stylepad(gPad);
        hRatio->Draw("COLZ");

        // Label the cut value at the top
        c->cd(0);
        TLatex lx;
        lx.SetNDC(); lx.SetTextSize(0.030);
        lx.DrawLatex(0.38, 0.965, Form("probRICH > %.2g", cut));

        c->Print(pdfName);

        // ── Save histograms in a per-cut directory ────────────────────────────
        outFile->mkdir(tag.c_str());
        outFile->cd(tag.c_str());
        hAll  ->Write();
        hPi   ->Write();
        hRatio->Write();

        delete hAll; delete hPi; delete hRatio;

        std::cout << "Cut " << cut << " done.\n";
    }

    c->Print(pdfName + "]");
    outFile->Close();
    std::cout << "PDF written to  : " << pdfName << "\n";
    std::cout << "ROOT written to : " << outName << "\n";
}
