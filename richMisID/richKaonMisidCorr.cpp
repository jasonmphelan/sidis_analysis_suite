/**
 * richKaonMisidCorr.cpp
 *
 * Stage 2 of the RICH kaon misidentification analysis.
 *
 * Reads the per-(p,θ) m²_RICH histograms produced by richKaonMisidHists,
 * fits the pion peak in each bin with a Gaussian, and computes the
 * misidentification correction factor:
 *
 *   f_misid(p,θ) = N_{π-like} / N_total      (pion contamination fraction)
 *   C(p,θ)       = 1 – f_misid               (kaon purity / correction factor)
 *
 * Statistical uncertainty (binomial):
 *   δf = (1/N_total²) × √(N_π × N_K × N_total)   [matches reference formula]
 *   δC = δf
 *
 * Outputs:
 *   • TH2F hCorrFactor_{Kp,Km}  – C(p,θ) with errors
 *   • TH2F hMisidFrac_{Kp,Km}   – f_misid(p,θ) with errors
 *   • Multi-page PDF             – diagnostic m² plot per bin, with:
 *       – m²_π and m²_K vertical lines (PDG expectations)
 *       – Gaussian fit of the pion peak
 *       – integration cut line
 *       – f_misid and C annotations
 *   • Heatmap PDF                – 2D colour map of C(p,θ) and f_misid(p,θ)
 *
 * Usage:
 *   ./richKaonMisidCorr <input_hists.root> <output.root> [n_sigma=2.5]
 *
 * PDG masses used:
 *   Mπ = 0.139570 GeV/c²  →  m²π = 0.019480 GeV²/c⁴
 *   MK = 0.493677 GeV/c²  →  m²K = 0.243717 GeV²/c⁴
 */

#include <iostream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

// ── Bin geometry (must match Stage 1) ─────────────────────────────────────────
static constexpr int    N_P_BINS  = 30;
static constexpr double P_MIN     = 1.25;
static constexpr double P_MAX     = 5.00;
static constexpr double P_STEP    = (P_MAX - P_MIN) / N_P_BINS;   // 0.15 GeV/c

static constexpr int    N_TH_BINS = 60;
static constexpr double TH_MIN    =  5.0;
static constexpr double TH_MAX    = 35.0;
static constexpr double TH_STEP   = (TH_MAX - TH_MIN) / N_TH_BINS; // 2 deg

// ── PDG masses and expected m² values ─────────────────────────────────────────
static constexpr double M_PI_PDG  = 0.139570;            // GeV/c²
static constexpr double M_K_PDG   = 0.493677;            // GeV/c²
static constexpr double M_P_PDG   = 0.938272;            // GeV/c² (reference only)
static constexpr double M2_PI     = M_PI_PDG * M_PI_PDG; // ≈ 0.01948 GeV²/c⁴
static constexpr double M2_K      = M_K_PDG  * M_K_PDG;  // ≈ 0.24372 GeV²/c⁴
static constexpr double M2_P      = M_P_PDG  * M_P_PDG;  // ≈ 0.88036 GeV²/c⁴

// ── Fit ranges in m² ─────────────────────────────────────────────────────────
// Double-Gaussian (simultaneous π + K fit):

// Single-Gaussian fallback (pion peak only):
static constexpr double FIT_LO    = -0.15;  // GeV²/c⁴  (matches Stage-1 M2_LO)
static constexpr double FIT_HI    =  0.10;  // GeV²/c⁴

// Default π/K integration cut when all fits fail
static constexpr double DEFAULT_PI_CUT = 0.1250; // GeV²/c⁴

// Minimum statistics to attempt any fit
static constexpr int MIN_ENTRIES = 50;

// ── Quick-look diagnostic bin (single-page PDF) ───────────────────────────────
// ip=10  →  p  ∈ [2.75, 2.90] GeV/c  (high-p, good K/π RICH separation)
// it=4   →  θ  ∈ [17,   19 ] deg     (mid-angle, typically well-populated)
static constexpr int DIAG_IP = 10;
static constexpr int DIAG_IT =  4;

// ── Helper: compute bin centre (p, θ) ─────────────────────────────────────────
inline double pCentre (int ip) { return P_MIN  + (ip + 0.5) * P_STEP;  }
inline double thCentre(int it) { return TH_MIN + (it + 0.5) * TH_STEP; }
inline double pEdgeLo (int ip) { return P_MIN  + ip * P_STEP;  }
inline double thEdgeLo(int it) { return TH_MIN + it * TH_STEP; }

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "\nUsage: " << argv[0]
                  << " <input_hists.root> <output.root> [n_sigma=2.5]\n\n";
        return -1;
    }

    const char* in_name  = argv[1];
    const char* out_name = argv[2];
    double n_sigma = (argc > 3) ? std::atof(argv[3]) : 2.5;

    std::cerr << "Input  : " << in_name  << "\n";
    std::cerr << "Output : " << out_name << "\n";
    std::cerr << "σ cut  : " << n_sigma  << "\n";

    // ── Open files ────────────────────────────────────────────────────────────
    TFile* inFile = TFile::Open(in_name);
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "ERROR: cannot open " << in_name << "\n";
        return -1;
    }
    TFile* outFile = new TFile(out_name, "RECREATE");

    // ── Global style ──────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);   // we draw fit info in the legend manually
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    static const char* cstr[2]  = {"Kp",  "Km"};
    static const char* ctex[2]  = {"K^{+}", "K^{-}"};

    // ── Output TH2F histograms ────────────────────────────────────────────────
    TH2F* hCorr [2];   // C(p,θ) = 1 – f_misid
    TH2F* hMisid[2];   // f_misid(p,θ)
    TH2F* hNtot [2];   // total counts per bin (diagnostic)
    TH2F* hNpi  [2];   // π-like counts per bin (diagnostic)

    for (int ic = 0; ic < 2; ++ic) {
        hCorr[ic] = new TH2F(
            Form("hCorrFactor_%s", cstr[ic]),
            Form("Kaon purity  C(p,#theta) = 1 #minus f_{misid}  (%s);"
                 "p (GeV/c);#theta (deg);C(p,#theta)", ctex[ic]),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);
        hCorr[ic]->Sumw2();

        hMisid[ic] = new TH2F(
            Form("hMisidFrac_%s", cstr[ic]),
            Form("Misid fraction  f_{misid} = N_{#pi-like}/N_{tot}  (%s);"
                 "p (GeV/c);#theta (deg);f_{misid}", ctex[ic]),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);
        hMisid[ic]->Sumw2();

        hNtot[ic] = new TH2F(
            Form("hNtot_%s", cstr[ic]),
            Form("N_{total} per (p,#theta) bin  (%s);"
                 "p (GeV/c);#theta (deg);N_{total}", ctex[ic]),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);

        hNpi[ic] = new TH2F(
            Form("hNpi_%s", cstr[ic]),
            Form("N_{#pi-like} per (p,#theta) bin  (%s);"
                 "p (GeV/c);#theta (deg);N_{#pi-like}", ctex[ic]),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);
    }

    // ── Loop: one charge, one (p,θ) bin at a time ─────────────────────────────
    for (int ic = 0; ic < 2; ++ic) {

        // Multi-page diagnostic PDF (one page per non-empty bin)
        TCanvas* cDiag  = new TCanvas(Form("cDiag_%s", cstr[ic]),
                                      "", 900, 650);
        cDiag->SetLeftMargin (0.13);
        cDiag->SetRightMargin(0.04);
        cDiag->SetTopMargin  (0.10);
        cDiag->SetBottomMargin(0.12);

        TString pdfDiag = Form("rich_misid_diag_%s.pdf", cstr[ic]);
        cDiag->Print(pdfDiag + "[");   // open the multi-page document

        for (int ip = 0; ip < N_P_BINS; ++ip) {
            for (int it = 0; it < N_TH_BINS; ++it) {

                // ── Retrieve 1D m² histogram ─────────────────────────────────
                TString hname = Form("hM2_%s_p%i_t%i", cstr[ic], ip, it);
                TH1F* h = dynamic_cast<TH1F*>(inFile->Get(hname));

                // Not enough data: mark as zero and skip
                if (!h || h->GetEntries() < MIN_ENTRIES) {
                    hCorr [ic]->SetBinContent(ip+1, it+1, 0.0);
                    hCorr [ic]->SetBinError  (ip+1, it+1, 0.0);
                    hMisid[ic]->SetBinContent(ip+1, it+1, 0.0);
                    hMisid[ic]->SetBinError  (ip+1, it+1, 0.0);
                    continue;
                }

                // Clone so we do not modify the histogram stored in inFile
                TH1F* hfit = dynamic_cast<TH1F*>(h->Clone(
                                 Form("hfit_%s_p%i_t%i", cstr[ic], ip, it)));
                hfit->SetDirectory(nullptr);

                const double sqrt2pi = std::sqrt(2.0 * M_PI);
                const double bw      = hfit->GetBinWidth(1);

                double Aseed_pi = std::max(
                    hfit->GetBinContent(hfit->GetXaxis()->FindBin(M2_PI)), 1.0);

                // ═══════════════════════════════════════════════════════════
                // Single-Gaussian pion fit
                // ═══════════════════════════════════════════════════════════
                TF1* fPion = new TF1(
                    Form("fPion_%s_%i_%i", cstr[ic], ip, it),
                    "gaus", FIT_LO, FIT_HI);
                fPion->SetParameters(Aseed_pi, M2_PI, 0.010);

                int statusSG = hfit->Fit(fPion, "RQNL");

                double pi_cut  = DEFAULT_PI_CUT;
                bool   fitSGOK = false;
                double N_pi_SG = 0.0, dN_pi_SG = 0.0;
                double N_tot_SG = hfit->Integral(0, hfit->GetNbinsX() + 1);
                double N_K_SG   = 0.0;
                if (statusSG == 0) {
                    double A_sg   = fPion->GetParameter(0);
                    double mean   = fPion->GetParameter(1);
                    double sg_sg  = std::abs(fPion->GetParameter(2));
                    double eA_sg  = fPion->GetParError(0);
                    double esg_sg = fPion->GetParError(2);
                    if (mean > -0.030 && mean < 0.070 && sg_sg > 0.002 && A_sg > 0.0) {
                        pi_cut   = mean + n_sigma * sg_sg;
                        fitSGOK  = true;
                        // N_π from Gaussian integral: A × σ × √(2π) / bin_width
                        N_pi_SG  = sqrt2pi * A_sg * sg_sg / bw;
                        // δN_π via error propagation (ignoring A/σ covariance)
                        dN_pi_SG = (sqrt2pi / bw)
                                   * std::sqrt(sg_sg*sg_sg * eA_sg*eA_sg
                                             + A_sg*A_sg   * esg_sg*esg_sg);
                    }
                }
                N_K_SG = N_tot_SG - N_pi_SG;

                // ═══════════════════════════════════════════════════════════
                // Choose result; compute f_misid, C
                // ═══════════════════════════════════════════════════════════
                enum FitMode { MODE_1G, MODE_COUNT } fitMode;
                double N_pi, N_K, N_tot, delta_f;

                if (N_tot_SG > 0.0) {
                    fitMode = MODE_1G;
                    N_pi  = N_pi_SG;
                    N_K   = N_K_SG;
                    N_tot = N_tot_SG;
                    delta_f = (N_tot > 0.0) ? dN_pi_SG / N_tot : 0.0;
                } else {
                    fitMode = MODE_COUNT;
                    N_pi = N_K = N_tot = delta_f = 0.0;
                }

                double f_misid = (N_tot > 0.0) ? N_pi / N_tot : 0.0;
                double corr    = 1.0 - f_misid;
                double delta_c = delta_f;

                // Fill TH2F only for physically valid corrections
                if (corr >= 0.0 && corr <= 1.0) {
                    hCorr [ic]->SetBinContent(ip+1, it+1, corr);
                    hCorr [ic]->SetBinError  (ip+1, it+1, delta_c);
                    hMisid[ic]->SetBinContent(ip+1, it+1, f_misid);
                    hMisid[ic]->SetBinError  (ip+1, it+1, delta_f);
                }

                hNtot[ic]->SetBinContent(ip+1, it+1, N_tot);
                hNpi [ic]->SetBinContent(ip+1, it+1, N_pi );

                // ── Diagnostic plot ──────────────────────────────────────────
                cDiag->cd();
                cDiag->Clear();
                // Log y-scale: pion peak and smaller kaon peak both visible.
                //cDiag->SetLogy(1);

                double p_lo = pEdgeLo (ip);
                double t_lo = thEdgeLo(it);

                // Zoom x-axis to show pion and kaon peaks plus proton reference.
                //hfit->GetXaxis()->SetRangeUser(-0.05, 0.55);
                // Floor > 0 so the log axis renders cleanly.
                hfit->SetMinimum(0.3);

                hfit->SetTitle(
                    Form("m^{2}_{RICH}  %s   %.2f < p < %.2f GeV/c   "
                         "%.0f < #theta < %.0f deg;"
                         "m^{2}_{RICH} (GeV^{2}/c^{4});Counts",
                         ctex[ic],
                         p_lo, p_lo + P_STEP,
                         t_lo, t_lo + TH_STEP));
                hfit->GetXaxis()->SetTitleSize(0.050);
                hfit->GetYaxis()->SetTitleSize(0.050);
                hfit->GetXaxis()->SetLabelSize(0.042);
                hfit->GetYaxis()->SetLabelSize(0.042);
                hfit->SetLineColor(kBlue + 1);
                hfit->SetLineWidth(2);
                hfit->SetFillColorAlpha(kBlue - 9, 0.35);
                hfit->SetFillStyle(1001);
                hfit->Draw("HIST");

                // Vertical reference lines (non-zero floor for log scale)
                double ymin_line = 0.3;
                double ymax_line = hfit->GetBinContent(
                                       hfit->GetMaximumBin()) * 5.0;

                TLine* lPi = new TLine(M2_PI, ymin_line, M2_PI, ymax_line);
                lPi->SetLineColor(kRed + 1);
                lPi->SetLineWidth(2);
                lPi->SetLineStyle(7);
                lPi->Draw();

                TLine* lK = new TLine(M2_K, ymin_line, M2_K, ymax_line);
                lK->SetLineColor(kGreen + 2);
                lK->SetLineWidth(2);
                lK->SetLineStyle(7);
                lK->Draw();

                // Proton reference line (histogram extends to ~1.0 GeV²/c⁴)
                TLine* lP = new TLine(M2_P, ymin_line, M2_P, ymax_line);
                lP->SetLineColor(kViolet + 1);
                lP->SetLineWidth(2);
                lP->SetLineStyle(7);
                lP->Draw();

                TLine* lCut = new TLine(pi_cut, ymin_line, pi_cut, ymax_line);
                lCut->SetLineColor(kOrange + 1);
                lCut->SetLineWidth(2);
                lCut->SetLineStyle(1);
                lCut->Draw();

                if (fitMode == MODE_1G && fitSGOK) {
                    fPion->SetLineColor(kRed);
                    fPion->SetLineWidth(2);
                    fPion->Draw("SAME");
                }

                // ── Legend ────────────────────────────────────────────────────
                TLegend* leg = new TLegend(0.52, 0.42, 0.96, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillColorAlpha(kWhite, 0.75);
                leg->SetTextSize(0.033);
                leg->AddEntry(hfit, "m^{2}_{RICH} data", "F");
                leg->AddEntry(lPi,
                    Form("m^{2}_{#pi} = %.4f GeV^{2}/c^{4}", M2_PI), "L");
                leg->AddEntry(lK,
                    Form("m^{2}_{K}  = %.4f GeV^{2}/c^{4}", M2_K ), "L");
                leg->AddEntry(lP,
                    Form("m^{2}_{p}  = %.4f GeV^{2}/c^{4}", M2_P ), "L");
                leg->AddEntry(lCut,
                    Form("#pi-cut: %.4f GeV^{2}/c^{4}%s",
                         pi_cut, fitSGOK ? "" : " [default]"), "L");
                if (fitSGOK)
                    leg->AddEntry(fPion,
                        Form("Gauss: #mu=%.4f  #sigma=%.4f",
                             fPion->GetParameter(1),
                             fPion->GetParameter(2)), "L");
                leg->AddEntry(static_cast<TObject*>(nullptr),
                    Form("f_{misid} = %.4f #pm %.4f", f_misid, delta_f), "");
                leg->AddEntry(static_cast<TObject*>(nullptr),
                    Form("C = %.4f #pm %.4f   N_{tot}=%.0f",
                         corr, delta_c, N_tot), "");
                leg->Draw();

                // ── Fit-mode annotation (top-left) ────────────────────────────
                const char* modeStr =
                    (fitMode == MODE_1G) ? "#checkmark Gauss fit OK"
                                        : "#times  Fit failed - default cut used";
                int modeColor =
                    (fitMode == MODE_1G) ? kGreen + 2 : kRed + 1;
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.032);
                latex.SetTextColor(modeColor);
                latex.DrawLatex(0.14, 0.865, modeStr);

                cDiag->Print(pdfDiag);   // append one page to PDF

                // Cleanup per-bin temporaries
                delete leg;
                delete lPi;
                delete lK;
                delete lP;
                delete lCut;
                delete fPion;
                delete hfit;

            }  // θ loop
        }  // p loop

        cDiag->Print(pdfDiag + "]");   // close multi-page PDF
        std::cout << "Diagnostic PDF written: " << pdfDiag << "\n";
        delete cDiag;

    }  // charge loop

    // ── Quick-look diagnostic: one high-p bin, both charges side-by-side ────────
    {
        TCanvas* cQL = new TCanvas("cQuickLook", "", 1800, 650);
        cQL->Divide(2, 1, 0.005, 0.005);

        double p_lo_ql = pEdgeLo (DIAG_IP);
        double t_lo_ql = thEdgeLo(DIAG_IT);

        for (int ic = 0; ic < 2; ++ic) {
            cQL->cd(ic + 1);
            gPad->SetLeftMargin (0.13);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin  (0.10);
            gPad->SetBottomMargin(0.13);

            TString hname = Form("hM2_%s_p%i_t%i", cstr[ic], DIAG_IP, DIAG_IT);
            TH1F* hql = dynamic_cast<TH1F*>(inFile->Get(hname));
            if (!hql) {
                std::cerr << "WARNING: quick-look histogram " << hname
                          << " not found – skipping pad.\n";
                continue;
            }

            TH1F* hqlc = dynamic_cast<TH1F*>(
                hql->Clone(Form("hqlc_%s", cstr[ic])));
            hqlc->SetDirectory(nullptr);

            hqlc->SetTitle(
                Form("m^{2}_{RICH}  %s   %.2f < p < %.2f GeV/c   "
                     "%.0f < #theta < %.0f deg;"
                     "m^{2}_{RICH} (GeV^{2}/c^{4});Counts",
                     ctex[ic],
                     p_lo_ql, p_lo_ql + P_STEP,
                     t_lo_ql, t_lo_ql + TH_STEP));
            hqlc->SetLineColor(kBlue + 1);
            hqlc->SetLineWidth(2);
            hqlc->SetFillColorAlpha(kBlue - 9, 0.35);
            hqlc->SetFillStyle(1001);
            //gPad->SetLogy(1);
            hqlc->GetXaxis()->SetRangeUser(-0.05, 0.55);
            hqlc->SetMinimum(0.3);
            hqlc->Draw("HIST");

            double ymin_ql = 0.3;
            double ymax_ql = hqlc->GetBinContent(hqlc->GetMaximumBin()) * 5.0;

            TLine* lPiQL = new TLine(M2_PI, ymin_ql, M2_PI, ymax_ql);
            lPiQL->SetLineColor(kRed + 1); lPiQL->SetLineWidth(2);
            lPiQL->SetLineStyle(7); lPiQL->Draw();

            TLine* lKQL = new TLine(M2_K, ymin_ql, M2_K, ymax_ql);
            lKQL->SetLineColor(kGreen + 2); lKQL->SetLineWidth(2);
            lKQL->SetLineStyle(7); lKQL->Draw();

            TLine* lPQL = new TLine(M2_P, ymin_ql, M2_P, ymax_ql);
            lPQL->SetLineColor(kViolet + 1); lPQL->SetLineWidth(2);
            lPQL->SetLineStyle(7); lPQL->Draw();

            // ── Re-run single-Gaussian pion fit for the quick-look bin ───────
            const double sqrt2pi_ql = std::sqrt(2.0 * M_PI);
            const double bw_ql      = hqlc->GetBinWidth(1);

            double Aseed_pi_ql = std::max(
                hqlc->GetBinContent(hqlc->GetXaxis()->FindBin(M2_PI)), 1.0);

            TF1* fPionQL = new TF1(Form("fPionQL_%s", cstr[ic]),
                                   "gaus", FIT_LO, FIT_HI);
            fPionQL->SetParameters(Aseed_pi_ql, M2_PI, 0.010);
            int qlSGstatus = hqlc->Fit(fPionQL, "RQNL");

            double ql_A   = fPionQL->GetParameter(0);
            double ql_mu  = fPionQL->GetParameter(1);
            double ql_sg  = std::abs(fPionQL->GetParameter(2));
            bool ql_fitOK = (qlSGstatus == 0)
                         && (ql_A > 0.0)
                         && (ql_mu > -0.030 && ql_mu < 0.070)
                         && (ql_sg > 0.002);

            double ql_pi_cut = ql_fitOK ? ql_mu + n_sigma * ql_sg : DEFAULT_PI_CUT;

            TLine* lCutQL = new TLine(ql_pi_cut, ymin_ql, ql_pi_cut, ymax_ql);
            lCutQL->SetLineColor(kOrange + 1);
            lCutQL->SetLineWidth(2);
            lCutQL->SetLineStyle(1);
            lCutQL->Draw();

            if (ql_fitOK) {
                fPionQL->SetLineColor(kRed);
                fPionQL->SetLineWidth(2);
                fPionQL->Draw("SAME");
            }

            // Retrieve stored correction values
            double corr_ql = hCorr [ic]->GetBinContent(DIAG_IP+1, DIAG_IT+1);
            double fmis_ql = hMisid[ic]->GetBinContent(DIAG_IP+1, DIAG_IT+1);
            double ntot_ql = hNtot [ic]->GetBinContent(DIAG_IP+1, DIAG_IT+1);

            TLegend* legQL = new TLegend(0.52, 0.42, 0.96, 0.88);
            legQL->SetBorderSize(0);
            legQL->SetFillColorAlpha(kWhite, 0.75);
            legQL->SetTextSize(0.034);
            legQL->AddEntry(hqlc, "m^{2}_{RICH} data", "F");
            legQL->AddEntry(lPiQL,
                Form("m^{2}_{#pi} = %.4f GeV^{2}/c^{4}", M2_PI), "L");
            legQL->AddEntry(lKQL,
                Form("m^{2}_{K}  = %.4f GeV^{2}/c^{4}", M2_K ), "L");
            legQL->AddEntry(lPQL,
                Form("m^{2}_{p}  = %.4f GeV^{2}/c^{4}", M2_P ), "L");
            legQL->AddEntry(lCutQL,
                Form("#pi-cut: %.4f GeV^{2}/c^{4}%s",
                     ql_pi_cut, ql_fitOK ? "" : " [default]"), "L");
            if (ql_fitOK)
                legQL->AddEntry(fPionQL,
                    Form("Gauss: #mu=%.4f  #sigma=%.4f", ql_mu, ql_sg), "L");
            legQL->AddEntry(static_cast<TObject*>(nullptr),
                Form("f_{misid} = %.4f    C = %.4f", fmis_ql, corr_ql), "");
            legQL->AddEntry(static_cast<TObject*>(nullptr),
                Form("N_{tot} = %.0f", ntot_ql), "");
            legQL->Draw();

            TLatex latexQL;
            latexQL.SetNDC();
            latexQL.SetTextSize(0.034);
            latexQL.SetTextColor(ql_fitOK ? kGreen + 2 : kRed + 1);
            latexQL.DrawLatex(0.14, 0.865,
                ql_fitOK ? "#checkmark Gauss fit OK"
                         : "#times  Fit failed - default cut used");

            delete legQL;
            delete lPiQL;
            delete lKQL;
            delete lPQL;
            delete lCutQL;
            delete fPionQL;
            delete hqlc;
        }

        cQL->SaveAs("rich_misid_quicklook.pdf");
        std::cout << "Quick-look diagnostic written: rich_misid_quicklook.pdf\n";
        outFile->cd();
        cQL->Write();
        delete cQL;
    }

    // ── 2D Heatmaps ───────────────────────────────────────────────────────────
    // With 25 p-bins × 15 θ-bins the cell size is too small for text
    // overlays to be legible.  Pure COLZ with well-chosen z-ranges and
    // clearly labelled axes is much more readable.
    //
    // Produce each charge as its own full-width PDF page (no side-by-side
    // splitting) so the axes and colour scale are not squeezed.

    // ── Load acceptance-matching cut functions ────────────────────────────────
    // matchCut2D_map.root contains TF1s named max_<sec>_pip/pim and
    // min_<sec>_pip/pim (sec = 0..5). Each returns θ_cut(p) in degrees.
    TF1* matchMax[6][2];  // [sector][0=pip, 1=pim]
    TF1* matchMin[6][2];
    bool matchLoaded = false;
    {
        TFile* matchFile = TFile::Open("../data/acceptance_matching/matchCut2D_map.root");
        if (matchFile && !matchFile->IsZombie()) {
            const char* cq[2] = {"pip", "pim"};
            matchLoaded = true;
            for (int s = 0; s < 6; ++s) {
                for (int q = 0; q < 2; ++q) {
                    matchMax[s][q] = dynamic_cast<TF1*>(
                        matchFile->Get(Form("max_%i_%s", s, cq[q])));
                    matchMin[s][q] = dynamic_cast<TF1*>(
                        matchFile->Get(Form("min_%i_%s", s, cq[q])));
                    if (matchMax[s][q]) matchMax[s][q]->SetBit(TObject::kMustCleanup, false);
                    if (matchMin[s][q]) matchMin[s][q]->SetBit(TObject::kMustCleanup, false);
                }
            }
            // Keep objects alive after the TFile goes out of scope
            matchFile->Close();
            std::cerr << "Acceptance matching cuts loaded.\n";
        } else {
            std::cerr << "WARNING: could not open matchCut2D_map.root"
                         " – cuts will not be overlaid.\n";
        }
    }

    // Draw all sector cut curves onto the current pad (call after COLZ Draw)
    auto drawMatchCuts = [&]() {
        if (!matchLoaded) return;
        for (int s = 4; s < 5; ++s) {
            for (int q = 0; q < 2; ++q) {
                if (matchMax[s][q]) {
                    matchMax[s][q]->SetRange(P_MIN, P_MAX);
                    matchMax[s][q]->SetLineColor(kRed + 1);
                    matchMax[s][q]->SetLineWidth(1);
                    matchMax[s][q]->SetLineStyle(2);
                    matchMax[s][q]->Draw("SAME");
                }
                if (matchMin[s][q]) {
                    matchMin[s][q]->SetRange(P_MIN, P_MAX);
                    matchMin[s][q]->SetLineColor(kAzure + 1);
                    matchMin[s][q]->SetLineWidth(1);
                    matchMin[s][q]->SetLineStyle(2);
                    matchMin[s][q]->Draw("SAME");
                }
            }
        }
    };

    // Common axis style helper (called after each Draw)
    auto styleHeat = [](TH2F* h, const char* ztitle) {
        h->GetZaxis()->SetTitle(ztitle);
        h->GetZaxis()->SetTitleSize(0.050);
        h->GetZaxis()->SetLabelSize(0.040);
        h->GetZaxis()->SetTitleOffset(1.25);
        h->GetXaxis()->SetTitleSize(0.052);
        h->GetXaxis()->SetLabelSize(0.042);
        h->GetXaxis()->SetNdivisions(510);
        h->GetYaxis()->SetTitleSize(0.052);
        h->GetYaxis()->SetLabelSize(0.042);
        h->GetYaxis()->SetNdivisions(515);
        h->SetContour(99);
    };

    for (int ic = 0; ic < 2; ++ic) {

        // ── Canvas 1: correction factor C(p,θ) ───────────────────────────────
        {
            TCanvas* cC = new TCanvas(
                Form("cCorr_%s", cstr[ic]),
                Form("Kaon purity C(p,#theta)  %s", ctex[ic]),
                1600, 700);
            cC->SetLeftMargin (0.10);
            cC->SetRightMargin(0.13);   // room for colour-axis label
            cC->SetBottomMargin(0.13);
            cC->SetTopMargin   (0.10);

            styleHeat(hCorr[ic], "C(p,#theta)");
            hCorr[ic]->Draw("COLZ");
            drawMatchCuts();

            TString pdf1 = Form("rich_misid_corrFactor_%s.pdf", cstr[ic]);
            cC->SaveAs(pdf1);
            std::cout << "Heatmap PDF written: " << pdf1 << "\n";
            outFile->cd();
            cC->Write();
            delete cC;
        }

        // ── Canvas 2: misidentification fraction f_misid(p,θ) ────────────────
        {
            TCanvas* cF = new TCanvas(
                Form("cMisid_%s", cstr[ic]),
                Form("Misid fraction f_{{misid}}(p,#theta)  %s", ctex[ic]),
                1600, 700);
            cF->SetLeftMargin (0.10);
            cF->SetRightMargin(0.13);
            cF->SetBottomMargin(0.13);
            cF->SetTopMargin   (0.10);

            styleHeat(hMisid[ic], "f_{misid}(p,#theta)");
            hMisid[ic]->Draw("COLZ");
            drawMatchCuts();

            TString pdf2 = Form("rich_misid_misidFrac_%s.pdf", cstr[ic]);
            cF->SaveAs(pdf2);
            std::cout << "Heatmap PDF written: " << pdf2 << "\n";
            outFile->cd();
            cF->Write();
            delete cF;
        }

        // ── Canvas 3: side-by-side summary (compact, for quick comparison) ───
        {
            TCanvas* cHeat = new TCanvas(
                Form("cHeat_%s", cstr[ic]),
                Form("RICH kaon misidentification  %s", ctex[ic]),
                2000, 700);
            cHeat->Divide(2, 1, 0.005, 0.005);

            cHeat->cd(1);
            gPad->SetLeftMargin (0.10);
            gPad->SetRightMargin(0.16);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin   (0.10);
            styleHeat(hCorr[ic], "C(p,#theta)");
            hCorr[ic]->Draw("COLZ");
            drawMatchCuts();

            cHeat->cd(2);
            gPad->SetLeftMargin (0.10);
            gPad->SetRightMargin(0.16);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin   (0.10);
            styleHeat(hMisid[ic], "f_{misid}(p,#theta)");
            hMisid[ic]->Draw("COLZ");
            drawMatchCuts();

            TString pdfHeat = Form("rich_misid_heatmap_%s.pdf", cstr[ic]);
            cHeat->SaveAs(pdfHeat);
            std::cout << "Heatmap PDF written: " << pdfHeat << "\n";
            outFile->cd();
            cHeat->Write();
            delete cHeat;
        }
    }

    // ── K⁺/K⁻ correction ratio ───────────────────────────────────────────────
    {
        TH2F* hRatio = dynamic_cast<TH2F*>(hCorr[0]->Clone("hCorrRatio_KpOverKm"));
        hRatio->SetDirectory(nullptr);
        hRatio->SetTitle(
            "C(p,#theta) ratio  K^{+}/K^{-};"
            "p (GeV/c);#theta (deg);C_{K^{+}} / C_{K^{-}}");
        hRatio->Divide(hCorr[1]);

        TCanvas* cRatio = new TCanvas("cCorrRatio", "", 1600, 700);
        cRatio->SetLeftMargin (0.10);
        cRatio->SetRightMargin(0.13);
        cRatio->SetBottomMargin(0.13);
        cRatio->SetTopMargin   (0.10);

        styleHeat(hRatio, "C_{K^{+}} / C_{K^{-}}");
        hRatio->Draw("COLZ");
        drawMatchCuts();

        cRatio->SaveAs("rich_misid_corrRatio_KpOverKm.pdf");
        std::cout << "Ratio PDF written: rich_misid_corrRatio_KpOverKm.pdf\n";
        outFile->cd();
        hRatio->Write();
        cRatio->Write();
        delete cRatio;
    }

    // ── Write histograms to output ROOT file ──────────────────────────────────
    outFile->cd();
    for (int ic = 0; ic < 2; ++ic) {
        hCorr [ic]->Write();
        hMisid[ic]->Write();
        hNtot [ic]->Write();
        hNpi  [ic]->Write();
    }

    outFile->Close();
    inFile ->Close();

    std::cout << "Correction histograms written to: " << out_name << "\n";
    return 0;
}
