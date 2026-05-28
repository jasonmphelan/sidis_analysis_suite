/**
 * richKaonMisidHists.cpp
 *
 * Stage 1 of the RICH kaon misidentification analysis.
 *
 * Reads the ePi TTree from kaon-tagged SIDIS skim files and fills
 * TH1F histograms of m²_RICH in a (p, θ) grid.  Output is written
 * to a ROOT file that is consumed by richKaonMisidCorr (Stage 2).
 *
 * m²_RICH = p² × (1 - β²_rich) / β²_rich
 *
 * Binning:
 *   p  : 25 bins, 1.25 – 5.00 GeV/c  (Δp = 0.15 GeV/c)
 *   θ  : 15 bins, 5.0  – 35.0 deg    (Δθ = 2 deg)
 *
 * Usage:
 *   ./richKaonMisidHists <output.root> [target: 0=RGB(default), 1=RGA]
 *                        [tree_path_override]
 *
 * The tree_path_override argument lets you point at a specific skim
 * file instead of the default set of kaon skims.
 *
 * Do NOT edit files in sidis_analysis_suite/kaonAnalysis.
 * This new analysis lives in sidis_analysis_suite/richMisid/.
 */

#include <iostream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

// Project headers (same include path as kaonAnalysis programs)
#include "pion.h"
#include "electron.h"
#include "constants.h"
#include "cut_values.h"

// ── Bin geometry (must match Stage 2) ────────────────────────────────────────
static constexpr int    N_P_BINS  = 30;
static constexpr double P_MIN     = 1.25;   // GeV/c
static constexpr double P_MAX     = 5.00;   // GeV/c
static constexpr double P_STEP    = (P_MAX - P_MIN) / N_P_BINS;   // 0.15 GeV/c

static constexpr int    N_TH_BINS = 60;
static constexpr double TH_MIN    =  5.0;   // degrees
static constexpr double TH_MAX    = 35.0;   // degrees
static constexpr double TH_STEP   = (TH_MAX - TH_MIN) / N_TH_BINS; // 2 deg

// ── m² histogram parameters ───────────────────────────────────────────────────
static constexpr int    N_M2_BINS = 215;
static constexpr double M2_LO     = -0.15;  // GeV²/c⁴
static constexpr double M2_HI     =  1.0;  // GeV²/c⁴  (covers π & K peaks)

using namespace constants;
using namespace cutVals;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "\nUsage: " << argv[0]
                  << " <output.root> [target: 0=RGB(default) 1=RGA]"
                     " [tree_path_override]\n\n";
        return -1;
    }

    const char* out_name    = argv[1];
    int         target      = (argc > 2) ? std::atoi(argv[2]) : 0;
    int         matching    = (argc > 3 )?std::atoi(argv[3]) : 0;
    const char* path_override = (argc > 3) ? argv[4] : nullptr;

    std::cerr << "Output : " << out_name << "\n";
    std::cerr << "Target : " << (target == 0 ? "RGB/deuterium" : "RGA/proton") << "\n";

    // ── Declare output histograms ────────────────────────────────────────────
    // Indexed by [charge] where 0 = K+, 1 = K–
    static const char* cstr[2] = {"Kp", "Km"};

    // Per-(p,θ) 1D m² distributions
    TH1F* hM2[2][N_P_BINS][N_TH_BINS];

    // 2D overview histograms (diagnostic)
    TH2F* hM2_2d[2];        // m²_RICH vs p (summed over θ)
    TH2F* hM2_theta_2d[2];  // m²_RICH vs θ (summed over p)
    TH2F* hM2TOF_2d[2];     // m²_TOF  vs p
    TH2F* hM2_vs_TOF[2];    // m²_RICH vs m²_TOF
    TH2F* hBeta_2d[2];      // β_RICH vs p
    TH2F* hTheta_p[2];      // θ vs p (acceptance footprint)
    TH2F* hPhi_p[2];        // φ vs p (azimuthal uniformity)
    TH1F* hSector[2];       // sector occupancy
    TH2F* hQ2_xB;           // Q² vs xB (electron kinematics)

    hQ2_xB = new TH2F(
        "hQ2_xB",
        "Q^{2} vs x_{B};x_{B};Q^{2} (GeV^{2})",
        50, 0.0, 0.8, 50, 1.0, 7.0);
    hQ2_xB->Sumw2();

    for (int ic = 0; ic < 2; ++ic) {
        hM2_2d[ic] = new TH2F(
            Form("hM2_p_%s", cstr[ic]),
            Form("m^{2}_{RICH} vs p  (%s);p (GeV/c);m^{2}_{RICH} (GeV^{2}/c^{4})", cstr[ic]),
            N_P_BINS, P_MIN, P_MAX, N_M2_BINS, M2_LO, M2_HI);
        hM2_2d[ic]->Sumw2();

        hM2_theta_2d[ic] = new TH2F(
            Form("hM2_theta_%s", cstr[ic]),
            Form("m^{2}_{RICH} vs #theta  (%s);#theta (deg);m^{2}_{RICH} (GeV^{2}/c^{4})", cstr[ic]),
            N_TH_BINS, TH_MIN, TH_MAX, N_M2_BINS, M2_LO, M2_HI);
        hM2_theta_2d[ic]->Sumw2();

        hM2TOF_2d[ic] = new TH2F(
            Form("hM2TOF_p_%s", cstr[ic]),
            Form("m^{2}_{TOF} vs p  (%s);p (GeV/c);m^{2}_{TOF} (GeV^{2}/c^{4})", cstr[ic]),
            N_P_BINS, P_MIN, P_MAX, N_M2_BINS, M2_LO, M2_HI);
        hM2TOF_2d[ic]->Sumw2();

        hM2_vs_TOF[ic] = new TH2F(
            Form("hM2_RICH_vs_TOF_%s", cstr[ic]),
            Form("m^{2}_{RICH} vs m^{2}_{TOF}  (%s);m^{2}_{TOF} (GeV^{2}/c^{4});m^{2}_{RICH} (GeV^{2}/c^{4})", cstr[ic]),
            N_M2_BINS, M2_LO, M2_HI, N_M2_BINS, M2_LO, M2_HI);
        hM2_vs_TOF[ic]->Sumw2();

        hBeta_2d[ic] = new TH2F(
            Form("hBeta_p_%s", cstr[ic]),
            Form("#beta_{RICH} vs p  (%s);p (GeV/c);#beta_{RICH}", cstr[ic]),
            N_P_BINS, P_MIN, P_MAX, 200, 0.85, 1.02);
        hBeta_2d[ic]->Sumw2();

        hTheta_p[ic] = new TH2F(
            Form("hTheta_p_%s", cstr[ic]),
            Form("#theta vs p  (%s, RICH tracks);p (GeV/c);#theta (deg)", cstr[ic]),
            N_P_BINS, P_MIN, P_MAX, N_TH_BINS, TH_MIN, TH_MAX);

        hPhi_p[ic] = new TH2F(
            Form("hPhi_p_%s", cstr[ic]),
            Form("#phi vs p  (%s, RICH tracks);p (GeV/c);#phi (deg)", cstr[ic]),
            N_P_BINS, P_MIN, P_MAX, 72, -180., 180.);

        hSector[ic] = new TH1F(
            Form("hSector_%s", cstr[ic]),
            Form("Sector occupancy  (%s, RICH tracks);Sector;Counts", cstr[ic]),
            6, 0.5, 6.5);

        for (int ip = 0; ip < N_P_BINS; ++ip) {
            for (int it = 0; it < N_TH_BINS; ++it) {
                double p_lo = P_MIN + ip * P_STEP;
                double p_hi = p_lo  + P_STEP;
                double t_lo = TH_MIN + it * TH_STEP;
                double t_hi = t_lo  + TH_STEP;

                hM2[ic][ip][it] = new TH1F(
                    Form("hM2_%s_p%i_t%i", cstr[ic], ip, it),
                    Form("m^{2}_{RICH}  (%s)  %.2f<p<%.2f GeV  %.0f<#theta<%.0f deg;"
                         "m^{2}_{RICH} (GeV^{2}/c^{4});Counts",
                         cstr[ic], p_lo, p_hi, t_lo, t_hi),
                    N_M2_BINS, M2_LO, M2_HI);
                hM2[ic][ip][it]->Sumw2();
            }
        }
    }

    // ── Build input TChain ───────────────────────────────────────────────────
    TChain* chain = new TChain("ePi");

    if (path_override) {
        std::cerr << "Using override path: " << path_override << "\n";
        chain->Add(path_override);
    } else {
        // Default kaon skim paths (relative to the build/richMisid directory,
        // same convention as the existing kaonAnalysis binaries)
        if (target == 0) {
            // RGB / deuterium – three beam energies
            chain->Add("../trees/final_skims/10.2/tight_skim_new.root");
            //chain->Add("../trees/final_skims/10.4/tight_skim.root");
            //chain->Add("../trees/final_skims/10.6/tight_skim.root");
        } else {
            // RGA / proton
            chain->Add("../trees/final_skims/10.2/final_skim_rga.root");
        }
    }

    long long totalEntries = chain->GetEntries();
    std::cerr << "TChain entries: " << totalEntries << "\n";
    if (totalEntries == 0) {
        std::cerr << "ERROR: no events found – check file paths.\n";
        return -1;
    }

    // ── TTreeReader ──────────────────────────────────────────────────────────
    TTreeReader reader(chain);

    // Electron branch (needed to satisfy the stored TTree layout even if
    // we do not use electron kinematics for this analysis)
    TTreeReaderValue<electron> e(reader, "e");

    // Pion / kaon-candidate branch
    TTreeReaderArray<pion> pi(reader, "pi");

    // Selection flags written by the skimmer
    TTreeReaderArray<bool> isGoodPion_no_acc(reader, "isGoodPion_no_acc");
    TTreeReaderArray<bool> isGoodPion(reader, "isGoodPion");


    // ── Event loop ───────────────────────────────────────────────────────────
    long long nProc = 0;
    long long nGoodTrack = 0;

    while (reader.Next()) {
        if (nProc % 200000 == 0)
            std::cout << "Events processed: " << nProc << "\n";
        ++nProc;

        int nPions = static_cast<int>(pi.end() - pi.begin());
        for (int i = 0; i < nPions; ++i) {

            // Basic quality selection (no acceptance-map cut; same as
            // makeKaonPlots with theta_cut == 0)
            if (matching==0  && !isGoodPion_no_acc[i]) continue;
            if (matching==1  && !isGoodPion[i]) continue;

            // Require a valid RICH hit
            double beta_rich = pi[i].getBeta_rich();
            if (beta_rich < 1.0e-4) continue;   // no RICH measurement
            //if (beta_rich >= 1.0)   continue;   // unphysical

            // Track kinematics
            double p_pi    = pi[i].get3Momentum().Mag();           // GeV/c
            double theta_pi = pi[i].get3Momentum().Theta() * rad_to_deg; // deg

            // Fiducial range checks
            if (p_pi    < P_MIN  || p_pi    >= P_MAX)   continue;
            if (theta_pi < TH_MIN || theta_pi >= TH_MAX) continue;

            // ── Compute m²_RICH ──────────────────────────────────────────────
            // m² = p² × (1 – β²) / β²  [cf. prompt formula using pi[i] API]
            double m2_RICH = p_pi * p_pi
                             * (1.0 - beta_rich * beta_rich)
                             / (beta_rich * beta_rich);

            // Charge index: 0 = positive (K+), 1 = negative (K–)
            int ic = (pi[i].getCharge() < 0) ? 1 : 0;

            // Bin indices
            int ip = static_cast<int>((p_pi    - P_MIN)  / P_STEP);
            int it = static_cast<int>((theta_pi - TH_MIN) / TH_STEP);
            // Clamp (should not be needed after range checks above)
            if (ip < 0 || ip >= N_P_BINS)   continue;
            if (it < 0 || it >= N_TH_BINS)  continue;

            // TOF m²
            double beta_tof = pi[i].getBeta();
            double m2_TOF = (beta_tof > 1.0e-4)
                ? p_pi * p_pi * (1.0 - beta_tof * beta_tof) / (beta_tof * beta_tof)
                : -999.;

            double phi_pi = pi[i].get3Momentum().Phi() * rad_to_deg;

            // Fill histograms
            hM2[ic][ip][it]->Fill(m2_RICH);
            hM2_2d[ic]->Fill(p_pi, m2_RICH);
            hM2_theta_2d[ic]->Fill(theta_pi, m2_RICH);
            if (beta_tof > 1.0e-4) {
                hM2TOF_2d[ic]->Fill(p_pi, m2_TOF);
                hM2_vs_TOF[ic]->Fill(m2_TOF, m2_RICH);
            }
            hBeta_2d[ic]->Fill(p_pi, beta_rich);
            hTheta_p[ic]->Fill(p_pi, theta_pi);
            hPhi_p[ic]->Fill(p_pi, phi_pi);
            hSector[ic]->Fill(pi[i].getSector());
            hQ2_xB->Fill(e->getXb(), e->getQ2());

            ++nGoodTrack;
        }
    }

    std::cout << "\nDone.  Events processed: " << nProc
              << "  RICH tracks accepted: " << nGoodTrack << "\n";

    // ── Write output ─────────────────────────────────────────────────────────
    TFile* outFile = new TFile(out_name, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Cannot create output file: " << out_name << "\n";
        return -1;
    }

    outFile->cd();
    hQ2_xB->Write();
    for (int ic = 0; ic < 2; ++ic) {
        hM2_2d[ic]->Write();
        hM2_theta_2d[ic]->Write();
        hM2TOF_2d[ic]->Write();
        hM2_vs_TOF[ic]->Write();
        hBeta_2d[ic]->Write();
        hTheta_p[ic]->Write();
        hPhi_p[ic]->Write();
        hSector[ic]->Write();
        for (int ip = 0; ip < N_P_BINS; ++ip) {
            for (int it = 0; it < N_TH_BINS; ++it) {
                // Only write bins with at least one entry to keep file lean
                if (hM2[ic][ip][it]->GetEntries() > 0)
                    hM2[ic][ip][it]->Write();
            }
        }
    }

    outFile->Close();
    std::cerr << "Histograms written to: " << out_name << "\n";

    delete chain;
    return 0;
}
