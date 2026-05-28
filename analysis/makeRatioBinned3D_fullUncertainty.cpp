// =============================================================================
// makeRatioBinned3D_fullUncertainty.cpp
//
// (a) Purpose
// -----------
// Drop-in companion to analysis/makeRatioBinned3D.cpp that recomputes the
// fragmentation-function ratio FF = (4 - R) / (4R - 1), R = sigma+/sigma-,
// from the same SIDIS event chains.  The central values for sigma+, sigma-,
// R, and FF are evaluated FULLY EVENT-BY-EVENT — every correction (C_mc,
// C_pi2k, C_k2pi, alpha_VM) is applied at each event's own kinematics
// during the fill loops, never at the bin centre.
//
// The statistical uncertainty uses the same scalar variance decomposition
// as the reference makeRatioBinned3D.cpp, generalised for event-by-event
// weights:
//
//   sigma^h = S_pi_evt^h - S_K_evt^h - S_W_evt^h - S_W_evt^{!h}
//
// where each S_X_evt^h is a per-event sum with C_mc folded in.  The
// variance keeps the reference's 5-term per-charge structure:
//
//   Var(sigma^h) = (B_raw^h)^2 * sigma_{C_mc^h}^2          (C_mc, same-charge)
//                + PoissonPID_pi^h + PoissonPID_K^h + Var_PID_C^h
//                                                          (counting + PID-corr)
//                + VMUnc^h                                  (VM, same-charge)
//                + (W_raw^{!h})^2 * sigma_{C_mc^{!h}}^2     (C_mc, cross-charge)
//                + VMUnc^{!h}                               (VM, cross-charge)
//
//   Cov(sigma+, sigma-) = -B_raw+ * W_raw+ * sigma_{C_mc+}^2
//                        -B_raw- * W_raw- * sigma_{C_mc-}^2
//                        +VMUnc+
//                        +VMUnc-
//
// Per-event generalisations of the reference quantities:
//   PoissonPID_pi^h = sum_i (C_mc,i^h * C_pi2k^h(p_i))^2       over pion events
//   PoissonPID_K^h  = sum_i (C_mc,i^h * C_k2pi^h(p_i))^2       over kaon events
//   Var_PID_C^h     = sum_p [ NMcPi^h(p)^2 * sigma_{C_pi2k(p)}^2
//                           + NMcK^h(p)^2  * sigma_{C_k2pi(p)}^2 ]
//     where NMcPi^h(p) = sum_{i in N_pi^h(p)} C_mc,i^h  (likewise NMcK^h)
//   VMUnc^h         = sum_i (C_mc,i^h * delta_alpha_i^h)^2     over VM events
//
//   B_raw^h = S_pi_raw^h - S_K_raw^h - S_W_raw^h
//   W_raw^h = S_W_raw^h
//     where S_*_raw^h are the analogous per-event sums WITHOUT C_mc
//     folded in (so they can be multiplied by sigma_{C_mc^h} evaluated at
//     the event-weighted bin-mean kinematics).
//
// (b) Reference
// -------------
// analysis/makeRatioBinned3D.cpp is the source of the structure, fill
// loops, central-value definitions, and the variance algebra ported here.
// THAT FILE IS NOT MODIFIED.  All helper symbols in this file are
// renamed with the suffix "_fullUnc" so the two compilation units can
// coexist without collisions.
//
// (c) Assumptions
// ---------------
// 1.  C_mc^h is treated as a single scalar with uncertainty sigma_{C_mc^h}
//     evaluated at the event-weighted bin-mean kinematics for the
//     same-/cross-charge C_mc term.  (Same convention as reference.)
// 2.  C_pi2k^h(p) and C_k2pi^h(p) are independent across p-bins and
//     across charges.  Their uncertainties enter Var_PID_C^h.
// 3.  Vector-meson weights alpha_i^h are independent across events and
//     charges; their per-event uncertainties delta_alpha_i^h are
//     uncorrelated and enter VMUnc^h.
// 4.  Raw counts N_pi^h(p), N_K^h(p) follow Poisson statistics; with
//     event-by-event weighting the counting variance is the sum of
//     squared weights, captured by PoissonPID_pi/K.
// 5.  Beam-energy bins are independent — variances add across E.
// =============================================================================

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <regex>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TRandom3.h"
#include "TParameter.h"
#include "electron.h"
#include "pion.h"
#include "constants.h"
#include "cut_values.h"
#include "correctionTools.h"
#include "analyzer.h"
#define CORR_PATH _DATA
#define HIST_PATH _HIST

using std::cerr;
using std::cout;
using std::isfinite;
using std::isnan;
using std::ofstream;
using namespace cutVals;
using namespace constants;

// -----------------------------------------------------------------------------
// Type aliases for the legacy (sample/charge/var/...) accumulators.
// Mirror analysis/makeRatioBinned3D.cpp so the fill loops are 1:1.
// -----------------------------------------------------------------------------
using VecE_fU      = std::vector<double>;
using VecP_fU      = std::vector<VecE_fU>;
using VecZ_fU      = std::vector<VecP_fU>;
using VecXB_fU     = std::vector<VecZ_fU>;
using VecQ2_fU     = std::vector<VecXB_fU>;
using VecVar_fU    = std::vector<VecQ2_fU>;
using VecCharge_fU = std::vector<VecVar_fU>;
using VecSample_fU = std::vector<VecCharge_fU>;

using VecZ_w_fU      = std::vector<VecE_fU>;
using VecXB_w_fU     = std::vector<VecZ_w_fU>;
using VecQ2_w_fU     = std::vector<VecXB_w_fU>;
using VecVar_w_fU    = std::vector<VecQ2_w_fU>;
using VecCharge_w_fU = std::vector<VecVar_w_fU>;
using VecSample_w_fU = std::vector<VecCharge_w_fU>;

// -----------------------------------------------------------------------------
// Type aliases for the scalar and per-p accumulators used by the
// full-event-by-event variance machinery.
//   ScalBin[charge][var][Q2][xB][Z][E]     — scalar per (charge,var,Q2,xB,Z,E)
//   PerPBin[charge][var][Q2][xB][Z][p][E]  — scalar per (...,p,E) for NMcPi/NMcK
// -----------------------------------------------------------------------------
using ScalE_fU       = std::vector<double>;            // [E]
using ScalZ_fU       = std::vector<ScalE_fU>;          // [Z][E]
using ScalXB_fU      = std::vector<ScalZ_fU>;          // [xB][Z][E]
using ScalQ2_fU      = std::vector<ScalXB_fU>;         // [Q2][xB][Z][E]
using ScalVar_fU     = std::vector<ScalQ2_fU>;         // [var][...]
using ScalChrg_fU    = std::vector<ScalVar_fU>;        // [charge][...]

using PerPE_fU       = std::vector<double>;            // [E]
using PerPP_fU       = std::vector<PerPE_fU>;          // [p][E]
using PerPZ_fU       = std::vector<PerPP_fU>;          // [Z][p][E]
using PerPXB_fU      = std::vector<PerPZ_fU>;          // [xB][...]
using PerPQ2_fU      = std::vector<PerPXB_fU>;         // [Q2][...]
using PerPVar_fU     = std::vector<PerPQ2_fU>;         // [var][...]
using PerPChrg_fU    = std::vector<PerPVar_fU>;        // [charge][...]

// -----------------------------------------------------------------------------
// Configuration globals (mirrors original; renamed-suffixed where needed).
// -----------------------------------------------------------------------------
int matchType_fU     = 0;
int applyCorr_fU     = 0;
int binnedOut_fU     = 1;
int map_fU           = 0;
int n4d_corr_bins_fU = 0;
int target_fU        = 0;
TString var_name_fU  = "null";
int bins_var_fU      = 1;
double var_min_fU    = 0;
double var_max_fU    = 1;
TString matching_file_fU = "matchCut2D_10.root";
int overflowBinning_fU = 0;  // bitmap: 1=Q2, 2=xB, 4=Z, 8=var
int verbose_fU       = 0;

// -----------------------------------------------------------------------------
// Forward declarations of helpers (renamed for symbol-coexistence).
// -----------------------------------------------------------------------------
double getVarVal_fullUnc(TString var, electron e, pion pi);
void   zeroSuppress_fullUnc(TH1F* h);
void   getNewRhoName_fullUnc(double E, TString& rho_norm_name);
bool   updateCorrectionsForBeam_fullUnc(double eBeam, double& beam_energy,
                                        TString& corrName, correctionTools& corr);
bool   updateCorrectionsForBeam_fullUnc(int E, TString& corrName,
                                        correctionTools& corr);

// -----------------------------------------------------------------------------
// calcMx2pi: missing mass of the X system recoiling against two pions.
//
//   Mx_2pi = | q + p_target - p_pi1 - p_pi2 |
//
// where q = virtual-photon 4-momentum = e.getQ(),
//       p_target = (0, 0, 0, Mp)  (target nucleon at rest),
//       p_pi1, p_pi2 = 4-momenta of the two pions.
//
// Returns Mx_2pi in GeV.  The caller should ensure both pions originate
// from the same event and apply any desired quality cuts before calling.
// -----------------------------------------------------------------------------
double calcMx2pi(const electron& e, const pion& pi1, const pion& pi2) {
    static const TLorentzVector p_target(0., 0., 0., Mp);
    TLorentzVector miss = e.getQ() + p_target - pi1.get4Momentum() - pi2.get4Momentum();
    return miss.Mag();
}

// Variance ingredient: PID-correction term (per charge).
//   Var_PID_C^h = sum_p [ NMcPi^h(p)^2 * sigma_{C_pi2k^h(p)}^2
//                       + NMcK^h(p)^2  * sigma_{C_k2pi^h(p)}^2 ]
static double computeVarPIDCorr_fullUnc(const std::vector<double>& NMcPi_p,
                                        const std::vector<double>& NMcK_p,
                                        const std::vector<double>& sigPi2k_p,
                                        const std::vector<double>& sigK2pi_p);

// FF + uncertainty packaging.
struct FFResult_fullUnc {
    double R, sigma_R;
    double FF, sigma_FF;
    bool   valid;          // false if sigma- == 0 / NaN, etc.
    bool   nearPole;       // true if |4R - 1| < 0.05
};
static FFResult_fullUnc computeFFAndUncertainty_fullUnc(
    double sp, double sm,
    double Var_sp, double Var_sm,
    double Cov_pm);

// Output filename helper: insert "_fullUnc" before the ".root" suffix.
static TString makeOutputName_fullUnc(TString in){
    if( in.EndsWith(".root") ){
        in.ReplaceAll(".root", "_fullUnc.root");
        return in;
    }
    return in + "_fullUnc.root";
}

// =============================================================================
// main()
// =============================================================================
int main( int argc, char** argv){

    if( argc < 13 ){
        cerr << "Usage: ./makeRatio3D_fullUnc [beam energy] [output file]\n";
        cerr << "       [match type] [apply corr] [map] [acceptance file]\n"
             << "       apply corr: bitwise — bits[1:0]: 0=none,1=bin migration,2=acceptance,3=full MC; bit2(4)=kaon k-factor; bit3(8)=rho/VM\n";
        cerr << "       [var name] [bins var] [var min] [var max] [binned out] [correction name]\n";
        cerr << "       (opt: Mx2 cut) (opt: rho mode) (opt: target) (opt: matching file)\n";
        cerr << "       (opt: overflow binning) (opt: verbose 0/1)\n";
        cerr << "       Output is written to <output file> with '_fullUnc' inserted before the .root extension.\n";
        return -1;
    }

    double inBeam = atof(argv[1]);

    TString out_name_raw = argv[2];
    matchType_fU = atoi(argv[3]);
    applyCorr_fU = atoi(argv[4]);
    map_fU       = atoi(argv[5]);
    TString acc_name = argv[6];

    var_name_fU = argv[7];
    bins_var_fU = atoi(argv[8]);
    var_min_fU  = atof(argv[9]);
    var_max_fU  = atof(argv[10]);
    if( atoi(argv[11]) > 1 ){ binnedOut_fU = bins_var_fU; }
    TString correction_name = argv[12];
    double Mx2_cut  = 1.25;
    double rho_mode = 1;
    if( argc > 13 ) Mx2_cut          = atof(argv[13]);
    if( argc > 14 ) rho_mode         = atoi(argv[14]);
    if( argc > 15 ) target_fU        = atoi(argv[15]);
    if( argc > 16 ) matching_file_fU = argv[16];
    if( argc > 17 ) overflowBinning_fU = atoi(argv[17]);
    if( argc > 18 ) verbose_fU       = atoi(argv[18]);

    TString out_name = makeOutputName_fullUnc(out_name_raw);

    cout << "Target: " << target_fU
         << " (" << (target_fU == 1 ? "RGA/proton" : "RGB/deuterium") << ")\n";
    if( overflowBinning_fU )
        cout << "Overflow binning bitmap=" << overflowBinning_fU
             << " (" << ((overflowBinning_fU & 1) ? "Q2 " : "")
             << ((overflowBinning_fU & 2) ? "xB " : "")
             << ((overflowBinning_fU & 4) ? "Z "  : "")
             << ((overflowBinning_fU & 8) ? "var" : "") << ")\n";

    // Probe correction file for 4D binning metadata before any array allocation
    if( bins_var_fU > 1 ){
        TString probePath = (TString)_DATA + "/correctionFiles/" + correction_name;
        TFile * pf = TFile::Open(probePath);
        if( pf && !pf->IsZombie() ){
            int count = 0;
            while( pf->Get(Form("hMcCorrP_fit_%d", count)) ) count++;
            n4d_corr_bins_fU = count;
            if( count > 0 ){
                auto* pMin = (TParameter<double>*)pf->Get("var_min");
                auto* pMax = (TParameter<double>*)pf->Get("var_max");
                if( pMin ) var_min_fU = pMin->GetVal();
                if( pMax ) var_max_fU = pMax->GetVal();
                cout << "Correction file: " << count << " 4D bins, var=["
                     << var_min_fU << ", " << var_max_fU << "]\n";
            } else {
                cout << "No 4D histograms in correction file; using 3D corrections (xB, Q2, z)\n";
            }
            delete pf;
        } else {
            cerr << "Warning: could not open correction file for probing: " << probePath << "\n";
        }
    }

    cout << "Beam energy: " << inBeam << "  Output: " << out_name << "\n";
    cout << "Match type: " << matchType_fU
         << "  Apply corr: " << applyCorr_fU << "  Map: " << map_fU << "\n";
    cout << "Var: " << var_name_fU
         << " [" << var_min_fU << ", " << var_max_fU << "] (" << bins_var_fU << " bins)\n";
    cout << "Corrections: " << correction_name
         << "  Mx2 cut: " << Mx2_cut << "\n";

    TFile * outFile = new TFile(out_name, "RECREATE");

    TChain * dChain = new TChain("ePi");
    TChain * kChain = new TChain("ePi");
    TChain * rChain = new TChain("ePi");

    TString base = "../trees/";
    constexpr int N_E = 3;

    if( target_fU == 0 ){
        if( inBeam == 0 || inBeam == 10.2){
            dChain->Add( base + "final_skims/10.2/final_skim.root");
            if( applyCorr_fU & 4 ) kChain->Add( base + "final_skims/kaons_10.2/final_skim.root");
            if( applyCorr_fU & 8 ) rChain->Add( base + "final_skims/rho_skims/rotated_10.2_final.root");
        }
        if( inBeam == 0 || inBeam == 10.4){
            dChain->Add( base + "final_skims/10.4/final_skim.root");
            if( applyCorr_fU & 4 ) kChain->Add( base + "final_skims/kaons_10.4/final_skim.root");
            if( applyCorr_fU & 8 ) rChain->Add( base + "final_skims/rho_skims/rotated_10.4_final.root");
        }
        if( inBeam == 0 || inBeam == 10.6){
            dChain->Add( base + "final_skims/10.6/final_skim.root");
            if( applyCorr_fU & 4 ) kChain->Add( base + "final_skims/kaons_10.6/final_skim.root");
            if( applyCorr_fU & 8 ) rChain->Add( base + "final_skims/rho_skims/rotated_10.6_final.root");
        }
    }
    if( target_fU == 1 ){
        dChain->Add( base + "final_skims/rga/final_skim_rga_comb.root");
        if( applyCorr_fU & 4 ) kChain->Add( base + "final_skims/kaons_10.2/final_skim_rga.root");
        if( applyCorr_fU & 8 ) rChain->Add( "../trees/rga/rotated_10.2_final.root");
    }

    cout << "Chains loaded\n";

    // -------------------------------
    // Dimensions
    // -------------------------------
    constexpr int nSamples = 4;
    constexpr int nCharge  = 2;
    int nVar = bins_var_fU;

    // -------------------------------
    // Legacy event/error/weight accumulators (preserved verbatim from original).
    // These drive the central-value calculation in the original fill loops.
    // -------------------------------
    VecSample_fU events_in_bin( nSamples, VecCharge_fU(
        nCharge, VecVar_fU( nVar, VecQ2_fU( bins_Q2,
            VecXB_fU( bins_xB, VecZ_fU( bins_Z,
                VecP_fU( bins_p, VecE_fU(N_E, 0.0) ) ) ) ) ) ) );

    VecSample_fU errors_in_bin( nSamples, VecCharge_fU(
        nCharge, VecVar_fU( nVar, VecQ2_fU( bins_Q2,
            VecXB_fU( bins_xB, VecZ_fU( bins_Z,
                VecP_fU( bins_p, VecE_fU(N_E, 0.0) ) ) ) ) ) ) );

    VecSample_w_fU weights_in_bin( nSamples, VecCharge_w_fU(
        nCharge, VecVar_w_fU( nVar, VecQ2_w_fU( bins_Q2,
            VecXB_w_fU( bins_xB, VecZ_w_fU( bins_Z,
                VecE_fU(N_E, 0.0) ) ) ) ) ) );

    // VM weight uncertainty sum (legacy, still used for fallback diagnostics).
    VecCharge_w_fU vm_weight_var_in_bin( nCharge, VecVar_w_fU(
        nVar, VecQ2_w_fU( bins_Q2, VecXB_w_fU( bins_xB,
            VecZ_w_fU( bins_Z, VecE_fU(N_E, 0.0) ) ) ) ) );

    // Event-weighted kinematic sums (used for mean correction lookup).
    VecVar_w_fU sum_Q2_kin( nVar, VecQ2_w_fU( bins_Q2, VecXB_w_fU( bins_xB,
        VecZ_w_fU( bins_Z, VecE_fU(N_E, 0.0) ) ) ) );
    VecVar_w_fU sum_xB_kin( nVar, VecQ2_w_fU( bins_Q2, VecXB_w_fU( bins_xB,
        VecZ_w_fU( bins_Z, VecE_fU(N_E, 0.0) ) ) ) );
    VecVar_w_fU sum_z_kin ( nVar, VecQ2_w_fU( bins_Q2, VecXB_w_fU( bins_xB,
        VecZ_w_fU( bins_Z, VecE_fU(N_E, 0.0) ) ) ) );
    VecVar_w_fU count_kin ( nVar, VecQ2_w_fU( bins_Q2, VecXB_w_fU( bins_xB,
        VecZ_w_fU( bins_Z, VecE_fU(N_E, 0.0) ) ) ) );

    // -------------------------------
    // correctionTools setup.
    // We MUST load the continuous corrections so that C_mc^h, C_pi2k^h(p),
    // C_k2pi^h(p), and their uncertainties are available per event from
    // the TF3/TF1 fits.
    // -------------------------------
    correctionTools corrector(n4d_corr_bins_fU > 0 ? 4 : 2);
    if( n4d_corr_bins_fU > 0 ){
        corrector.setN4dBins(n4d_corr_bins_fU);
    }
    corrector.setWeightFitName( correction_name );
    corrector.setMisidFitName ( "corrections_misid_fit.root" );
    corrector.loadContinuousCorrections();
    corrector.printFilePaths();

    // -------------------------------
    // dOu heatmaps: TH2D(Z, Q2) for each xB bin, indexed ix2..ix13.
    // hCombo_ix[ixB] covers x in [0.3,0.7] (Z) and y in [2,8] (Q2).
    // Bins are stored 0-indexed here: heatmap_dOu[i] corresponds to hCombo_ix(i+2).
    // -------------------------------
    constexpr int N_DOu_IX   = 12;   // ix2 .. ix13
    constexpr int DOu_IX_OFF =  2;   // first index in the file
    TH2D* heatmap_dOu[N_DOu_IX] = {};
    {
        TString hmPath = (TString)_DATA + "/corrections/dOu_heatmaps.root";
        TFile* hmFile = TFile::Open(hmPath);
        if( !hmFile || hmFile->IsZombie() ){
            cerr << "ERROR: could not open dOu heatmap file: " << hmPath << "\n";
        } else {
            for( int k = 0; k < N_DOu_IX; ++k ){
                TString name = Form("hCombo_ix%d", k + DOu_IX_OFF);
                heatmap_dOu[k] = (TH2D*)hmFile->Get(name);
                if( !heatmap_dOu[k] ){
                    cerr << "WARNING: missing heatmap " << name << " in " << hmPath << "\n";
                } else {
                    heatmap_dOu[k]->SetDirectory(nullptr);  // detach from file
                }
            }
            hmFile->Close();
            delete hmFile;
            cout << "Loaded " << N_DOu_IX << " dOu heatmaps from " << hmPath << "\n";
        }
    }

    // bits[1:0] of applyCorr_fU: 0=none → mc_factor_type=-1 (no MC corrections needed)
    const int mc_factor_type = (applyCorr_fU & 3) - 1;

    // -------------------------------
    // Scalar and per-p accumulators for the full-uncertainty machinery.
    // All sums are event-by-event with the appropriate per-event weights.
    //
    //   S_pi_evt[c][var][Q2][xB][Z][E]  = sum_{i in N_pi^c}  C_mc,i^c * C_pi2k^c(p_i)
    //   S_K_evt [c][...]                = sum_{i in N_K^c}   C_mc,i^c * C_k2pi^c(p_i)
    //   S_W_evt [c][...]                = sum_{i in n_c}     C_mc,i^c * alpha_i^c
    //   (the three above give the central sigma^h — see calc loop)
    //
    //   S_pi_raw[c][...]                = sum_{i in N_pi^c}  C_pi2k^c(p_i)
    //   S_K_raw [c][...]                = sum_{i in N_K^c}   C_k2pi^c(p_i)
    //   S_W_raw [c][...]                = sum_{i in n_c}     alpha_i^c
    //   (raw sums, no C_mc — used in the (B^h)^2 * sigma_{C_mc^h}^2 term
    //    and Cov terms that factor out C_mc^h to its bin-mean value.)
    //
    //   PoissonPID_pi[c][...]           = sum (C_mc,i^c * C_pi2k^c(p_i))^2
    //   PoissonPID_K [c][...]           = sum (C_mc,i^c * C_k2pi^c(p_i))^2
    //     (Poisson counting variance generalised by sum-of-squared-weights.)
    //
    //   VMUnc[c][...]                   = sum (C_mc,i^c * delta_alpha_i^c)^2
    //     (Per-event VM weight uncertainty with C_mc folded in.)
    //
    //   NMcPi[c][var][Q2][xB][Z][p][E]  = sum_{i in N_pi^c(p)} C_mc,i^c
    //   NMcK [c][var][Q2][xB][Z][p][E]  = sum_{i in N_K^c(p)}  C_mc,i^c
    //     (Used by Var_PID_C^h: per-p MC-weighted counts that multiply
    //      sigma_{C_pi2k^h(p)} and sigma_{C_k2pi^h(p)}.)
    // -------------------------------
    auto makeScal = [&](){
        return ScalChrg_fU( nCharge, ScalVar_fU( nVar,
            ScalQ2_fU( bins_Q2, ScalXB_fU( bins_xB,
                ScalZ_fU( bins_Z, ScalE_fU(N_E, 0.0) ) ) ) ) );
    };
    auto makePerP = [&](){
        return PerPChrg_fU( nCharge, PerPVar_fU( nVar,
            PerPQ2_fU( bins_Q2, PerPXB_fU( bins_xB,
                PerPZ_fU( bins_Z, PerPP_fU( bins_p,
                    PerPE_fU(N_E, 0.0) ) ) ) ) ) );
    };

    ScalChrg_fU S_pi_evt = makeScal();
    ScalChrg_fU S_K_evt  = makeScal();
    ScalChrg_fU S_W_evt  = makeScal();

    ScalChrg_fU S_pi_raw = makeScal();
    ScalChrg_fU S_K_raw  = makeScal();
    ScalChrg_fU S_W_raw  = makeScal();

    ScalChrg_fU PoissonPID_pi = makeScal();
    ScalChrg_fU PoissonPID_K  = makeScal();
    ScalChrg_fU VMUnc         = makeScal();

    PerPChrg_fU NMcPi = makePerP();
    PerPChrg_fU NMcK  = makePerP();

    // -------------------------------
    // Output histograms (mirrors original).
    // -------------------------------
    TH1F * hZ[bins_var_fU][bins_Q2][bins_xB][2];
    TH1F * hZ_k[bins_var_fU][bins_Q2][bins_xB][2];
    TH1F * hZ_r[bins_var_fU][bins_Q2][bins_xB][2];
    TString charge_str[2] = {"", "_Pim"};

    cout << "Make histograms\n";
    for (int var = 0; var < bins_var_fU; var++) {
        for (int i = 0; i < bins_Q2; i++) {
            for (int j = 0; j < bins_xB; j++) {
                for (int k = 0; k < 2; k++) {
                    TString suffix = (bins_var_fU > 1)
                        ? Form("_%d_%d_%d", i+1, j+1, k+1)
                        : Form("_%d_%d", i+1, j+1);
                    TString hbase = Form("_%d%s%s", var, charge_str[k].Data(), suffix.Data());
                    hZ  [var][i][j][k] = new TH1F("hRatio"  + hbase,
                                                  Form("hRatio %s%s", charge_str[k].Data(), suffix.Data()),
                                                  bins_Z, 0.3, 1.0);
                    hZ_k[var][i][j][k] = new TH1F("hRatio_k"+ hbase,
                                                  Form("hRatio_k %s%s", charge_str[k].Data(), suffix.Data()),
                                                  bins_Z, 0.3, 1.0);
                    hZ_r[var][i][j][k] = new TH1F("hRatio_r"+ hbase,
                                                  Form("hRatio_r %s%s", charge_str[k].Data(), suffix.Data()),
                                                  bins_Z, 0.3, 1.0);
                    hZ  [var][i][j][k]->Sumw2();
                    hZ_k[var][i][j][k]->Sumw2();
                    hZ_r[var][i][j][k]->Sumw2();
                }
            }
        }
    }

    analyzer anal( 0, -1 );
    anal.setAnalyzerLevel(0);
    anal.loadMatchingFunctions(matching_file_fU);
    anal.loadMatchingFunctions3D();
    anal.loadAcceptanceMapContinuous( (TString)_DATA + (TString)"/acceptance_map/" + acc_name );


    // -----------------------------------------------------------------
    // Pion analysis loop (mirror of original; with full-unc accumulators).
    // -----------------------------------------------------------------
    cout << "Begin pion analysis\n";
    double beam_energy = 10.2;

    TTreeReader reader_rec(dChain);
    TTreeReaderValue<double>      eBeam(reader_rec, "Ebeam");
    TTreeReaderValue<electron>    e    (reader_rec, "e");
    TTreeReaderArray<pion>        pi   (reader_rec, "pi");
    TTreeReaderArray<bool>        isGoodPion_no_acc(reader_rec, "isGoodPion_no_acc");
    TTreeReaderArray<bool>        isGoodPion       (reader_rec, "isGoodPion");
    TTreeReaderArray<bool>        isGoodPion3D     (reader_rec, "isGoodPion_3d");

    struct Bins_fU { int E, Q2, xB, Z, var, p; };
    auto calcBins = [&](const electron& el, const pion& pi_ref, double p_mom) -> Bins_fU {
        int bE   = (int)((beam_energy - 10.2) / .2);
        int bQ2  = (int)(((el.getQ2()  - Q2_min) / (Q2_max - Q2_min)) * bins_Q2);
        int bxB  = (int)(((el.getXb()  - xB_min) / (xB_max - xB_min)) * bins_xB);
        int bZ   = (int)(((pi_ref.getZ() - .3)    / (1. - .3))         * bins_Z);
        int bVar = (int)(((getVarVal_fullUnc(var_name_fU, el, pi_ref) - var_min_fU)
                          / (var_max_fU - var_min_fU)) * bins_var_fU);
        int bP   = -1;
        for (int bin = 0; bin < bins_p; bin++)
            if (p_mom > p_bin_edges[bin] && p_mom < p_bin_edges[bin+1]) { bP = bin; break; }
        if (bQ2  < 0 || (!(overflowBinning_fU & 1) && bQ2  >= bins_Q2))      bQ2  = -1;
        if (bxB  < 0 || (!(overflowBinning_fU & 2) && bxB  >= bins_xB))      bxB  = -1;
        if (bZ   < 0 || (!(overflowBinning_fU & 4) && bZ   >= bins_Z))        bZ   = -1;
        if (bVar < 0 || (!(overflowBinning_fU & 8) && bVar >= bins_var_fU))  bVar = -1;
        if (bE   < 0 || bE  >= N_E) bE = -1;
        if ((overflowBinning_fU & 1) && bQ2  >= bins_Q2)     bQ2  = bins_Q2     - 1;
        if ((overflowBinning_fU & 2) && bxB  >= bins_xB)     bxB  = bins_xB     - 1;
        if ((overflowBinning_fU & 4) && bZ   >= bins_Z)       bZ   = bins_Z      - 1;
        if ((overflowBinning_fU & 8) && bVar >= bins_var_fU) bVar = bins_var_fU - 1;
        return {bE, bQ2, bxB, bZ, bVar, bP};
    };

    auto applyMatching = [&](const pion& pi_ref, double p_mom) -> bool {
        int chargeIdx = (int)( pi_ref.getCharge() < 1 )+1;
        if (matchType_fU == 1 || matchType_fU == 2)
            return anal.applyAcceptanceMatching(pi_ref, 2);
        int chargeIdx_1 = ( matchType_fU == 3 ) ? 1 : chargeIdx;
        int chargeIdx_2 = ( matchType_fU == 3 ) ? 2 : chargeIdx;
        if ( matchType_fU == 3 || matchType_fU == 4 )
            return  anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 1, chargeIdx_1) >= 0 &&
                    anal.applyAcceptanceMap(p_mom, rad_to_deg*pi_ref.get3Momentum().Phi(), rad_to_deg*pi_ref.get3Momentum().Theta(), 2, chargeIdx_2) >= 0 && anal.applyAcceptanceMatching(pi_ref, 2);
        return true;
    };

    int event_total = reader_rec.GetEntries();
    while (reader_rec.Next()) {
        int event_count = reader_rec.GetCurrentEntry();
        if(event_count%100000 == 0)
            cout << "Pions: " << event_count << " / " << event_total << "\n";

        updateCorrectionsForBeam_fullUnc(*eBeam, beam_energy, correction_name, corrector);

        TVector3 e_mom = e->get3Momentum();
        if(matchType_fU<=1 && map_fU
           && anal.applyAcceptanceMap(e_mom.Mag(),
                                       rad_to_deg*e_mom.Phi(),
                                       rad_to_deg*e_mom.Theta(), 0) < 0 ) continue;

        for( int i = 0; i < (int)( pi.end() - pi.begin() ); i++ ){
            int chargeIdx = (int)( pi[i].getCharge() < 1 );
            int pid       = abs(pi[i].getPID_eb());
            if( pid > 321 ) continue;
            if( pi[i].getMx() < 1.7 ) continue;
            int pid_bin   = (pid==211) ? 0 : 1;
            double p_pi   = pi[i].get3Momentum().Mag();
            //if( (int)( pi.end() - pi.begin() ) == 2 && pow(calcMx2pi(*e, pi[0], pi[1]) , 2) < Mx2_cut) continue;
            if(matchType_fU==0 && !isGoodPion_no_acc[i]) continue;
            else if(matchType_fU==1 && !isGoodPion_no_acc[i]) continue;
            else if(matchType_fU==2 && !isGoodPion_no_acc[i]) continue;
            else if(matchType_fU==3 && (!isGoodPion_no_acc[i])) continue;
            else if(matchType_fU==4 && (!isGoodPion_no_acc[i])) continue;

            if(matchType_fU<=1 && map_fU
               && anal.applyAcceptanceMap( p_pi,
                                           rad_to_deg*pi[i].get3Momentum().Phi(),
                                           rad_to_deg*pi[i].get3Momentum().Theta(),
                                           chargeIdx + 1 ) < 0) continue;

            if( !applyMatching(pi[i], p_pi) ) continue;

            auto b = calcBins(*e, pi[i], p_pi);
            if (b.E < 0 || b.Q2 < 0 || b.xB < 0 || b.Z < 0 || b.var < 0) continue;

            // Per-event kinematics for the corrections.
            double xB_evt = e->getXb();
            double Q2_evt = e->getQ2();
            double z_evt  = pi[i].getZ();

            // dOu cut: require heatmap value in (0, 1).
            {
                int hm_idx = b.xB - DOu_IX_OFF;
                if (hm_idx >= 0 && hm_idx < N_DOu_IX && heatmap_dOu[hm_idx]) {
                    TH2D* h = heatmap_dOu[hm_idx];
                    double z_lk  = std::max(h->GetXaxis()->GetXmin() + 1e-9,
                                   std::min(z_evt,  h->GetXaxis()->GetXmax() - 1e-9));
                    double q2_lk = std::max(h->GetYaxis()->GetXmin() + 1e-9,
                                   std::min(Q2_evt, h->GetYaxis()->GetXmax() - 1e-9));
                    double dOu_val = h->Interpolate(z_lk, q2_lk);
                    if (dOu_val <= 0.0 || dOu_val >= 1.0) continue;
                }
            }

            corrector.setKinematics(xB_evt, Q2_evt, z_evt, p_pi);
            if( n4d_corr_bins_fU > 0 ) corrector.set4dBin(b.var);

            // Legacy weight (preserves original central-value structure).
            double weight = 1;
            if( (applyCorr_fU & 3) == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
            if( (applyCorr_fU & 3) == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
            if( (applyCorr_fU & 3) == 3 ) weight *= corrector.getCorrectionFactor(2, chargeIdx);
            if( (applyCorr_fU & 4) && b.p >= 0 )
                weight *= (pid==211)
                    ? corrector.getCorrectionFactor(3, chargeIdx)
                    : corrector.getCorrectionFactor(4, chargeIdx);

            // Per-event C_mc^h (used in the full-unc accumulators).
            double C_mc_evt = (mc_factor_type >= 0)
                ? corrector.getCorrectionFactor(mc_factor_type, chargeIdx) : 1.0;

            // PID factors per event (per p, per charge, per species).
            double C_pi2k = (applyCorr_fU & 4) ? corrector.getCorrectionFactor(3, chargeIdx) : 1.0;
            double C_k2pi = (applyCorr_fU & 4) ? corrector.getCorrectionFactor(4, chargeIdx) : 0.0;

            int Q2_lo  = (overflowBinning_fU & 1) ? 0 : b.Q2;
            int Q2_hi  = b.Q2 + 1;
            int xB_lo  = (overflowBinning_fU & 2) ? 0 : b.xB;
            int xB_hi  = b.xB + 1;
            int Z_lo   = (overflowBinning_fU & 4) ? 0 : b.Z;
            int Z_hi   = b.Z  + 1;
            int var_lo = (overflowBinning_fU & 8) ? 0 : b.var;
            int var_hi = b.var + 1;

            for (int iq = Q2_lo; iq < Q2_hi; iq++) {
                for (int ix = xB_lo; ix < xB_hi; ix++) {
                    for (int iz = Z_lo; iz < Z_hi; iz++) {
                        for (int ivar = var_lo; ivar < var_hi; ivar++) {
                        if( !(applyCorr_fU & 4) && pid != 211 ) continue;

                        // Legacy fills (preserves central-value structure).
                        events_in_bin [pid_bin][chargeIdx][ivar][iq][ix][iz][b.p][b.E]++;
                        weights_in_bin[pid_bin][chargeIdx][ivar][iq][ix][iz][b.E] += weight;
                        errors_in_bin [pid_bin][chargeIdx][ivar][iq][ix][iz][b.p][b.E]++;

                        sum_Q2_kin[ivar][iq][ix][iz][b.E] += Q2_evt;
                        sum_xB_kin[ivar][iq][ix][iz][b.E] += xB_evt;
                        sum_z_kin [ivar][iq][ix][iz][b.E] += z_evt;
                        count_kin [ivar][iq][ix][iz][b.E] += 1.0;

                        // === Full-unc accumulators (per-event sums) ===
                        if (b.p >= 0) {
                            if (pid == 211) {
                                // Pion-tagged event: contributes to S_PID via +C_pi2k * N_pi.
                                double w_pi = C_mc_evt * C_pi2k;
                                PoissonPID_pi[chargeIdx][ivar][iq][ix][iz][b.E] += w_pi * w_pi;
                                NMcPi[chargeIdx][ivar][iq][ix][iz][b.p][b.E]    += C_mc_evt;
                                // Event-by-event central-value sum (with C_mc folded in).
                                S_pi_evt[chargeIdx][ivar][iq][ix][iz][b.E]      += w_pi;
                                // Raw sum (no C_mc) used in the C_mc-sigma term of Var.
                                S_pi_raw[chargeIdx][ivar][iq][ix][iz][b.E]      += C_pi2k;
                            } else if (applyCorr_fU & 4) {
                                // RICH-tagged kaon: contributes to S_PID via -C_k2pi * N_K.
                                double w_k = C_mc_evt * C_k2pi;
                                PoissonPID_K[chargeIdx][ivar][iq][ix][iz][b.E] += w_k * w_k;
                                NMcK[chargeIdx][ivar][iq][ix][iz][b.p][b.E]    += C_mc_evt;
                                // Event-by-event central-value sum (with C_mc folded in).
                                S_K_evt[chargeIdx][ivar][iq][ix][iz][b.E]      += w_k;
                                // Raw sum (no C_mc).
                                S_K_raw[chargeIdx][ivar][iq][ix][iz][b.E]      += C_k2pi;
                            }
                        }
                        } // ivar
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // Rho / VM analysis loop (mirror of original; with full-unc fills).
    // -----------------------------------------------------------------
    TTreeReader reader_r(rChain);
    TTreeReaderValue<electron>      e_r(reader_r, "e");
    TTreeReaderArray<pion>          r  (reader_r, "pi");
    TTreeReaderArray<double>        rhoWeight(reader_r, "rhoWeight");
    TTreeReaderValue<double>        Mx_2pi   (reader_r, "Mx_2pi");
    TTreeReaderArray<double>        rhoError (reader_r, "rhoErr");
    TTreeReaderArray<bool>          isGoodRho(reader_r, "isGoodPion");
    TTreeReaderValue<TLorentzVector> beam    (reader_r, "beam");

    event_total = reader_r.GetEntries();

    if( applyCorr_fU & 8 ){
        cout << "Begin rho analysis\n";
        while (reader_r.Next()) {
            int event_count = reader_r.GetCurrentEntry();
            if(event_count%100000 == 0)
                cout << "Rhos:  " << event_count << " / " << event_total << "\n";

            updateCorrectionsForBeam_fullUnc(beam->E(), beam_energy, correction_name, corrector);

            int i = -1;
            if      ( isGoodRho[rho_mode] && r[rho_mode].getMx() > 1.7  ) i = rho_mode;
            else if ( rho_mode==1 && !isGoodRho[1] && r[0].getMx()>1.7
                      && r[1].get3Momentum().Theta()*rad_to_deg
                         > anal.getMaxTheta(r[1].get3Momentum().Mag(), 1) ) i = 0;
            else if ( rho_mode==0 && !isGoodRho[0] && r[1].getMx()>1.7
                      && r[0].get3Momentum().Theta()*rad_to_deg
                         < anal.getMinTheta(r[0].get3Momentum().Mag(), 2) ) i = 1;
            else { continue; }

            int chargeIdx = (int)( r[i].getCharge() < 1 );
            double p_pi   = r[i].get3Momentum().Mag();

            if( !applyMatching(r[i], p_pi) ) continue;
            if( (*Mx_2pi)*(*Mx_2pi) > Mx2_cut ) continue;

            const int index = 2;
            auto b = calcBins(*e_r, r[i], p_pi);
            if (b.E < 0 || b.Q2 < 0 || b.xB < 0 || b.Z < 0 || b.var < 0) continue;

            double xB_evt = e_r->getXb();
            double Q2_evt = e_r->getQ2();
            double z_evt  = r[i].getZ();

            // dOu cut: require heatmap value in (0, 1).
            {
                int hm_idx = b.xB - DOu_IX_OFF;
                if (hm_idx >= 0 && hm_idx < N_DOu_IX && heatmap_dOu[hm_idx]) {
                    TH2D* h = heatmap_dOu[hm_idx];
                    double z_lk  = std::max(h->GetXaxis()->GetXmin() + 1e-9,
                                   std::min(z_evt,  h->GetXaxis()->GetXmax() - 1e-9));
                    double q2_lk = std::max(h->GetYaxis()->GetXmin() + 1e-9,
                                   std::min(Q2_evt, h->GetYaxis()->GetXmax() - 1e-9));
                    double dOu_val = h->Interpolate(z_lk, q2_lk);
                    if (dOu_val <= 0.0 || dOu_val >= 1.0) continue;
                }
            }

            corrector.setKinematics(xB_evt, Q2_evt, z_evt, p_pi);
            if( n4d_corr_bins_fU > 0 ) corrector.set4dBin(b.var);

            double weight = rhoWeight[i];
            if( (applyCorr_fU & 3) == 1 ) weight *= corrector.getCorrectionFactor(0, chargeIdx);
            if( (applyCorr_fU & 3) == 2 ) weight *= corrector.getCorrectionFactor(1, chargeIdx);
            if( (applyCorr_fU & 3) == 3 ) weight *= corrector.getCorrectionFactor(2, chargeIdx);

            // Per-event C_mc^h for the full-unc accumulators.
            double C_mc_evt = (mc_factor_type >= 0)
                ? corrector.getCorrectionFactor(mc_factor_type, chargeIdx) : 1.0;

            double alpha_i  = rhoWeight[i];
            double dAlpha_i = rhoError [i];

            int Q2_lo  = (overflowBinning_fU & 1) ? 0 : b.Q2;
            int Q2_hi  = b.Q2 + 1;
            int xB_lo  = (overflowBinning_fU & 2) ? 0 : b.xB;
            int xB_hi  = b.xB + 1;
            int Z_lo   = (overflowBinning_fU & 4) ? 0 : b.Z;
            int Z_hi   = b.Z  + 1;
            int var_lo = (overflowBinning_fU & 8) ? 0 : b.var;
            int var_hi = b.var + 1;

            for (int iq = Q2_lo; iq < Q2_hi; iq++) {
                for (int ix = xB_lo; ix < xB_hi; ix++) {
                    for (int iz = Z_lo; iz < Z_hi; iz++) {
                        for (int ivar = var_lo; ivar < var_hi; ivar++) {
                        // Legacy fills (preserves central value structure).
                        events_in_bin [index][chargeIdx][ivar][iq][ix][iz][b.p][b.E] +=  rhoWeight[i];//(rhoWeight[i]-1);
                        weights_in_bin[index][chargeIdx][ivar][iq][ix][iz][b.E]      += weight;
                        errors_in_bin [index][chargeIdx][ivar][iq][ix][iz][b.p][b.E] += pow(rhoWeight[i],2) + pow(rhoWeight[i],2);//pow(rhoWeight[i]-1,2) + pow(rhoError[i]*(rhoWeight[i]-1),2);
                        vm_weight_var_in_bin[chargeIdx][ivar][iq][ix][iz][b.E]      += pow(rhoError[i],2);//pow(rhoError[i]*(rhoWeight[i]-1),2);

                        // === Full-unc accumulators (VM term) ===
                        double w_du = C_mc_evt * dAlpha_i;
                        VMUnc[chargeIdx][ivar][iq][ix][iz][b.E]    += w_du * w_du;
                        // Event-by-event central-value sum (with C_mc folded in).
                        S_W_evt[chargeIdx][ivar][iq][ix][iz][b.E]  += C_mc_evt * alpha_i;
                        // Raw sum (no C_mc).
                        S_W_raw[chargeIdx][ivar][iq][ix][iz][b.E]  += alpha_i;
                        } // ivar
                    }
                }
            }
        }
    }

    //////////////////////////////////////
    ///////////// Calcs //////////////////
    //////////////////////////////////////
    //
    // Full analytical statistical uncertainty on the fragmentation
    // function FF = (4 - R) / (4R - 1), where R = sigma+/sigma-, with
    // ALL central-value corrections applied event-by-event.
    //
    // Per-bin cross-section definition (per charge h in {+,-}):
    //
    //   sigma^h = S_pi_evt^h - S_K_evt^h - S_W_evt^h - S_W_evt^{!h}
    //
    // where each S_X_evt^h is a per-event sum (corrections folded in
    // at every event's actual kinematics during the fill loops above).
    //
    // The five-term variance decomposition (per charge) matches the
    // reference makeRatioBinned3D.cpp, with each ingredient generalised
    // for event-by-event weighting:
    //   [0] = (B_raw^h)^2 * sigma_{C_mc^h}^2        (C_mc, same-charge)
    //   [1] = PoissonPID_pi^h + PoissonPID_K^h
    //         + Var_PID_C^h                         (counting + PID-corr)
    //   [2] = VMUnc^h                                (VM, same-charge)
    //   [3] = (W_raw^{!h})^2 * sigma_{C_mc^{!h}}^2   (C_mc, cross-charge)
    //   [4] = VMUnc^{!h}                             (VM, cross-charge)
    //
    // Cov(sigma+, sigma-) four-term decomposition:
    //   [0] = -B_raw+ * W_raw+ * sigma_{C_mc+}^2
    //   [1] = -B_raw- * W_raw- * sigma_{C_mc-}^2
    //   [2] = +VMUnc+
    //   [3] = +VMUnc-
    //
    // Var(R) and Var(FF) then follow from standard error propagation
    // (formulas (6) and (7) in the reference).
    //

    cout << "Start calculations\n";

    // ── Diagnostic histograms for uncertainty decomposition ──
    // hSigma[var][Q2][xB][charge]: sigma+/- with delta-sigma as error bar
    // hCovHist[var][Q2][xB]:       Cov(sigma+, sigma-) per z-bin
    // hFF[var][Q2][xB]:            FF = (4 - R)/(4R - 1)
    TH1F * hSigma   [bins_var_fU][bins_Q2][bins_xB][2];
    TH1F * hCovHist [bins_var_fU][bins_Q2][bins_xB];
    TH1F * hFF      [bins_var_fU][bins_Q2][bins_xB];

    for (int var = 0; var < bins_var_fU; var++) {
        for (int ii = 0; ii < bins_Q2; ii++) {
            for (int jj = 0; jj < bins_xB; jj++) {
                TString ds = Form("_%d_%d_%d", var, ii+1, jj+1);
                for (int c = 0; c < 2; c++) {
                    hSigma[var][ii][jj][c] = new TH1F(
                        Form("hSigma%s%s%s", (c==0 ? "Plus" : "Minus"),
                             charge_str[c].Data(), ds.Data()),
                        "", bins_Z, 0.3, 1.0);
                    hSigma[var][ii][jj][c]->Sumw2();
                }
                hCovHist[var][ii][jj] = new TH1F(
                    Form("hCov%s", ds.Data()), "", bins_Z, 0.3, 1.0);
                hFF[var][ii][jj] = new TH1F(
                    Form("hFF%s", ds.Data()),
                    ";z;FF = (4 - R)/(4R - 1)",
                    bins_Z, 0.3, 1.0);
                hFF[var][ii][jj]->Sumw2();
            }
        }
    }

    double xB_step = (xB_max - xB_min) / (double)bins_xB;
    double Q2_step = (Q2_max - Q2_min) / (double)bins_Q2;

    // ── Determine which beam energies have data ──
    bool beamHasData[N_E] = {false, false, false};
    for (int E = 0; E < N_E; E++)
        for (int c = 0; c < 2 && !beamHasData[E]; c++)
            for (int v = 0; v < bins_var_fU && !beamHasData[E]; v++)
                for (int q = 0; q < bins_Q2 && !beamHasData[E]; q++)
                    for (int x = 0; x < bins_xB && !beamHasData[E]; x++)
                        for (int z = 0; z < bins_Z && !beamHasData[E]; z++)
                            for (int p = 0; p < bins_p && !beamHasData[E]; p++)
                                if (events_in_bin[0][c][v][q][x][z][p][E] > 0)
                                    beamHasData[E] = true;

    for (int E = 0; E < N_E; E++)
        cout << "Beam energy E=" << E << " (Ebeam=" << 10.2 + 0.2*E
             << "): " << (beamHasData[E] ? "HAS DATA" : "no data, skipping") << "\n";

    // ── Pre-load corrections once per beam energy ──
    correctionTools* corrByE[N_E] = {};
    for (int E = 0; E < N_E; E++) {
        if (!beamHasData[E]) continue;
        corrByE[E] = new correctionTools(n4d_corr_bins_fU > 0 ? 4 : 2);
        corrByE[E]->setN4dBins(n4d_corr_bins_fU);
        corrByE[E]->setMisidFitName("corrections_misid_fit.root");
        TString corrName_e = correction_name;
        updateCorrectionsForBeam_fullUnc(E, corrName_e, *corrByE[E]);
        cout << "Pre-loaded corrections for E=" << E
             << " (Ebeam=" << 10.2 + 0.2*E << ")\n";
    }

    // ── Main calculation loop ──
    for (int i = 1; i <= bins_Q2; i++) {
        for (int j = 1; j <= bins_xB; j++) {
            for (int var = 0; var < bins_var_fU; var++) {

                double xB_center = xB_min + (j - 0.5) * xB_step;
                double Q2_center = Q2_min + (i - 0.5) * Q2_step;

                // Per-z-bin accumulators (summed over beam energies E)
                double sigma_total[2][bins_Z];        // sigma+ and sigma-
                double dsigma_sq_total[2][bins_Z];    // (delta-sigma+-)^2
                double cov_total_z[bins_Z];           // Cov(sigma+, sigma-)

                // Per-term breakdown of Var(sigma^h) [c][t][z].  Indices:
                //   [0] = (B_raw^h)^2 * sigma_{C_mc^h}^2        (C_mc, same-charge)
                //   [1] = PoissonPID_pi^h + PoissonPID_K^h
                //         + Var_PID_C^h                         (counting + PID-corr)
                //   [2] = VMUnc^h                                (VM, same-charge)
                //   [3] = (W_raw^{!h})^2 * sigma_{C_mc^{!h}}^2   (C_mc, cross-charge)
                //   [4] = VMUnc^{!h}                             (VM, cross-charge)
                double dsig_term[2][5][bins_Z];

                // Four-source categorised breakdown of Var(sigma^h) [c][cat][z]:
                //   [0] counting-statistics : PoissonPID_pi^h + PoissonPID_K^h
                //   [1] PID-correction      : Var_PID_C^h
                //   [2] VM-weight           : VMUnc^h + VMUnc^{!h}
                //   [3] C_mc                : (B_raw^h)^2  * sigma_{C_mc^h}^2
                //                           + (W_raw^{!h})^2 * sigma_{C_mc^{!h}}^2
                double dsig_cat[2][4][bins_Z];

                // Covariance term breakdown [t][z]:
                //   [0] = -B_raw+ * W_raw+ * sigma_{C_mc+}^2
                //   [1] = -B_raw- * W_raw- * sigma_{C_mc-}^2
                //   [2] = +VMUnc+
                //   [3] = +VMUnc-
                double cov_term[4][bins_Z];

                memset(sigma_total,     0, sizeof(sigma_total));
                memset(dsigma_sq_total, 0, sizeof(dsigma_sq_total));
                memset(cov_total_z,     0, sizeof(cov_total_z));
                memset(dsig_term,       0, sizeof(dsig_term));
                memset(dsig_cat,        0, sizeof(dsig_cat));
                memset(cov_term,        0, sizeof(cov_term));

                for (int E = 0; E < N_E; E++) {
                    if (!beamHasData[E] || !corrByE[E]) continue;
                    correctionTools& corrE = *corrByE[E];

                    for (int zbin = 0; zbin < bins_Z; zbin++) {
                        double z_center = hZ[var][i-1][j-1][0]->GetXaxis()->GetBinCenter(zbin+1);

                        // Use event-weighted mean kinematics when overflow
                        // binning is active; fall back to bin centre otherwise.
                        double cnt   = count_kin[var][i-1][j-1][zbin][E];
                        double Q2_kin = (cnt > 0) ? sum_Q2_kin[var][i-1][j-1][zbin][E] / cnt : Q2_center;
                        double xB_kin = (cnt > 0) ? sum_xB_kin[var][i-1][j-1][zbin][E] / cnt : xB_center;
                        double z_kin  = (cnt > 0) ? sum_z_kin [var][i-1][j-1][zbin][E] / cnt : z_center;

                        // ── Get C_MC+- and delta-C_MC+- at the event-weighted bin-mean ──
                        // These multiply the "raw" (no-C_mc) sums in the C_mc-sigma terms,
                        // and replace the reference's bin-centre C_mc.  All other
                        // corrections (C_pi2k, C_k2pi, alpha_VM, plus C_mc itself) have
                        // already been folded in event-by-event during the fill loops.
                        corrE.setKinematics(xB_kin, Q2_kin, z_kin, 1.5);
                        if (n4d_corr_bins_fU > 0) corrE.set4dBin(var);

                        double C_MC[2], dC_MC[2];
                        for (int c = 0; c < 2; c++) {
                            C_MC[c]  = (mc_factor_type >= 0) ?
                                corrE.getCorrectionFactor(mc_factor_type, c) : 1.0;
                            dC_MC[c] = (mc_factor_type >= 0) ?
                                corrE.getCorrectionError(mc_factor_type, c) : 0.0;
                        }
                        (void)C_MC;  // C_MC itself is not used in the variance
                                      // (it's already in the event-by-event sums);
                                      // we only need its uncertainty dC_MC.

                        // =====================================================
                        // Central values (event-by-event, no bin-centre approximations):
                        //   sigma^h = S_pi_evt^h - S_K_evt^h - S_W_evt^h - S_W_evt^{!h}
                        // =====================================================
                        double Spi_eb[2], SK_eb[2], SW_eb[2];
                        for (int c = 0; c < 2; c++) {
                            Spi_eb[c] = S_pi_evt[c][var][i-1][j-1][zbin][E];
                            SK_eb [c] = S_K_evt [c][var][i-1][j-1][zbin][E];
                            SW_eb [c] = S_W_evt [c][var][i-1][j-1][zbin][E];
                        }
                        double sig[2];
                        sig[0] = Spi_eb[0] - SK_eb[0] - SW_eb[0] - SW_eb[1];
                        sig[1] = Spi_eb[1] - SK_eb[1] - SW_eb[1] - SW_eb[0];

                        // =====================================================
                        // "Raw" (no-C_mc) sums for the C_mc-sigma variance terms.
                        // These are the per-event analogues of B^h and W^h from the
                        // reference (which used bin-centre C_mc), with C_mc factored
                        // out so it can be multiplied by sigma_{C_mc^h} at the
                        // bin-mean kinematics.
                        //   B_raw^h = S_pi_raw^h - S_K_raw^h - S_W_raw^h
                        //   W_raw^h = S_W_raw^h
                        // =====================================================
                        double Spi_rw[2], SK_rw[2], SW_rw[2], B_raw[2], W_raw[2];
                        for (int c = 0; c < 2; c++) {
                            Spi_rw[c] = S_pi_raw[c][var][i-1][j-1][zbin][E];
                            SK_rw [c] = S_K_raw [c][var][i-1][j-1][zbin][E];
                            SW_rw [c] = S_W_raw [c][var][i-1][j-1][zbin][E];
                            B_raw [c] = Spi_rw[c] - SK_rw[c] - SW_rw[c];
                            W_raw [c] = SW_rw[c];
                        }

                        // =====================================================
                        // Var_PID_C^h: PID-correction uncertainty propagated through
                        //   sum_p NMcPi^h(p)^2 * sigma_{C_pi2k^h(p)}^2
                        //        + NMcK^h(p)^2  * sigma_{C_k2pi^h(p)}^2
                        // (The per-p PID sigmas are the only quantities still
                        //  evaluated at bin-mean kinematics — the misid TF1
                        //  uncertainty fit is intrinsically per-p.)
                        // =====================================================
                        std::vector<double> NMcPi_p[2], NMcK_p[2], sigPi2k_p[2], sigK2pi_p[2];
                        for (int c = 0; c < 2; c++) {
                            NMcPi_p[c].assign(bins_p, 0.0);
                            NMcK_p [c].assign(bins_p, 0.0);
                            sigPi2k_p[c].assign(bins_p, 0.0);
                            sigK2pi_p[c].assign(bins_p, 0.0);
                            for (int pbin = 0; pbin < bins_p; pbin++) {
                                NMcPi_p[c][pbin] = NMcPi[c][var][i-1][j-1][zbin][pbin][E];
                                NMcK_p [c][pbin] = NMcK [c][var][i-1][j-1][zbin][pbin][E];
                                double p_center = (p_bin_edges[pbin] + p_bin_edges[pbin+1]) / 2.;
                                corrE.setKinematics(xB_kin, Q2_kin, z_kin, p_center);
                                if (n4d_corr_bins_fU > 0) corrE.set4dBin(var);
                                sigPi2k_p[c][pbin] = (applyCorr_fU & 4)
                                    ? corrE.getCorrectionError(3, c) : 0.0;
                                sigK2pi_p[c][pbin] = (applyCorr_fU & 4)
                                    ? corrE.getCorrectionError(4, c) : 0.0;
                            }
                        }

                        // =====================================================
                        // Collect the remaining per-charge variance ingredients
                        // already pre-accumulated event-by-event during the fills.
                        // =====================================================
                        double PPID_pi_v[2], PPID_K_v[2], VarPIDC[2], VMU[2];
                        for (int c = 0; c < 2; c++) {
                            PPID_pi_v[c] = PoissonPID_pi[c][var][i-1][j-1][zbin][E];
                            PPID_K_v [c] = PoissonPID_K [c][var][i-1][j-1][zbin][E];
                            VarPIDC  [c] = computeVarPIDCorr_fullUnc(
                                NMcPi_p[c], NMcK_p[c], sigPi2k_p[c], sigK2pi_p[c]);
                            VMU      [c] = VMUnc[c][var][i-1][j-1][zbin][E];
                        }

                        // =====================================================
                        // Var(sigma+) and Var(sigma-) — reference's 5-term form
                        // generalised for event-by-event weighting.
                        //
                        //   t0[c][0] = (B_raw^h)^2 * sigma_{C_mc^h}^2     (C_mc same)
                        //   t0[c][1] = PoissonPID_pi^h + PoissonPID_K^h
                        //              + Var_PID_C^h                     (counting+PID)
                        //   t0[c][2] = VMUnc^h                            (VM same)
                        //   t0[c][3] = (W_raw^{!h})^2 * sigma_{C_mc^{!h}}^2 (C_mc cross)
                        //   t0[c][4] = VMUnc^{!h}                         (VM cross)
                        // =====================================================
                        double t0[2][5];
                        // sigma+ terms (c = 0)
                        t0[0][0] = pow(B_raw[0], 2) * pow(dC_MC[0], 2);
                        t0[0][1] = PPID_pi_v[0] + PPID_K_v[0] + VarPIDC[0];
                        t0[0][2] = VMU[0];
                        t0[0][3] = pow(W_raw[1], 2) * pow(dC_MC[1], 2);
                        t0[0][4] = VMU[1];
                        // sigma- terms (c = 1)
                        t0[1][0] = pow(B_raw[1], 2) * pow(dC_MC[1], 2);
                        t0[1][1] = PPID_pi_v[1] + PPID_K_v[1] + VarPIDC[1];
                        t0[1][2] = VMU[1];
                        t0[1][3] = pow(W_raw[0], 2) * pow(dC_MC[0], 2);
                        t0[1][4] = VMU[0];

                        double dsig_sq[2] = {0, 0};
                        for (int c = 0; c < 2; c++)
                            for (int t = 0; t < 5; t++) {
                                dsig_sq[c]            += t0[c][t];
                                dsig_term[c][t][zbin] += t0[c][t];
                            }

                        // Four-source categorised breakdown of Var(sigma^h).
                        for (int c = 0; c < 2; c++) {
                            int oc = 1 - c;  // opposite-charge index
                            dsig_cat[c][0][zbin] += PPID_pi_v[c] + PPID_K_v[c];
                            dsig_cat[c][1][zbin] += VarPIDC[c];
                            dsig_cat[c][2][zbin] += VMU[c] + VMU[oc];
                            dsig_cat[c][3][zbin] += pow(B_raw[c],  2) * pow(dC_MC[c],  2)
                                                  + pow(W_raw[oc], 2) * pow(dC_MC[oc], 2);
                        }

                        // =====================================================
                        // Cov(sigma+, sigma-) — reference's four-term form
                        // with per-event weighting.
                        // =====================================================
                        double ct[4];
                        ct[0] = -B_raw[0] * W_raw[0] * pow(dC_MC[0], 2);
                        ct[1] = -B_raw[1] * W_raw[1] * pow(dC_MC[1], 2);
                        ct[2] =  VMU[0];
                        ct[3] =  VMU[1];

                        double cov_E = ct[0] + ct[1] + ct[2] + ct[3];

                        // ── Accumulate across beam energies (independent E => Var adds) ──
                        for (int c = 0; c < 2; c++) {
                            sigma_total    [c][zbin] += sig[c];
                            dsigma_sq_total[c][zbin] += dsig_sq[c];
                        }
                        cov_total_z[zbin] += cov_E;
                        for (int t = 0; t < 4; t++)
                            cov_term[t][zbin] += ct[t];

                    } // zbin
                } // E

                // ══════════════════════════════════════════════════════════════
                // SET BIN CONTENT — output histograms filled here
                //
                // Inputs come from the E-loop above (sigma_total, dsigma_sq_total,
                // cov_total_z), which accumulated contributions from each beam energy.
                // ══════════════════════════════════════════════════════════════
                for (int zbin = 0; zbin < bins_Z; zbin++) {
                    double sp     = sigma_total[0][zbin];
                    double sm     = sigma_total[1][zbin];
                    double dsp_sq = dsigma_sq_total[0][zbin];
                    double dsm_sq = dsigma_sq_total[1][zbin];
                    double cov    = cov_total_z[zbin];

                    // ── Bin content: sigma+/- (cross-sections, with statistical error) ──
                    hSigma[var][i-1][j-1][0]->SetBinContent(zbin+1, sp);
                    hSigma[var][i-1][j-1][0]->SetBinError  (zbin+1, sqrt(fmax(dsp_sq, 0.0)));
                    hSigma[var][i-1][j-1][1]->SetBinContent(zbin+1, sm);
                    hSigma[var][i-1][j-1][1]->SetBinError  (zbin+1, sqrt(fmax(dsm_sq, 0.0)));
                    // ── Bin content: Cov(sigma+, sigma-) ──
                    hCovHist[var][i-1][j-1]->SetBinContent(zbin+1, cov);

                    // ── R = sigma+/sigma- and FF = (4 - R)/(4R - 1) with full propagation ──
                    auto res = computeFFAndUncertainty_fullUnc(sp, sm, dsp_sq, dsm_sq, cov);
                    if (res.valid) {
                        hZ[var][i-1][j-1][0]->SetBinContent(zbin+1, res.R);
                        hZ[var][i-1][j-1][0]->SetBinError  (zbin+1, res.sigma_R);
                        hFF[var][i-1][j-1]->SetBinContent(zbin+1, res.FF);
                        hFF[var][i-1][j-1]->SetBinError  (zbin+1, res.sigma_FF);
                        if (res.nearPole) {
                            cerr << "WARNING: |4R-1| < 0.05 at "
                                 << "(Q2=" << i << ", xB=" << j
                                 << ", var=" << var << ", z=" << zbin
                                 << "), R = " << res.R
                                 << " — Taylor expansion for sigma_FF may be unreliable.\n";
                        }
                    } else {
                        hZ[var][i-1][j-1][0]->SetBinContent(zbin+1, 0);
                        hZ[var][i-1][j-1][0]->SetBinError  (zbin+1, 0);
                        hFF[var][i-1][j-1]->SetBinContent(zbin+1, 0);
                        hFF[var][i-1][j-1]->SetBinError  (zbin+1, 0);
                    }
                } // zbin

                // ── Diagnostic: fractional contribution of each variance source ──
                // Gated behind the `verbose` flag (CLI arg 18).  Prints, per bin,
                // the fractional contribution to Var(sigma+) and Var(sigma-) of
                // the four sources: counting statistics, PID correction,
                // VM weights, and C_mc; plus the covariance and (dR/R)^2
                // breakdowns for completeness.
                if (verbose_fU) {
                    // Accumulate term sums across z-bins
                    double sum_dsig[2]     = {0, 0};
                    double sum_term[2][5]  = {{0}};
                    double sum_cat[2][4]   = {{0}};
                    double sum_ct[4]       = {0};
                    double sum_dRoR_sq     = 0;
                    double sum_dRoR_sp = 0, sum_dRoR_sm = 0, sum_dRoR_cov = 0;
                    int nz_filled = 0;

                    for (int zbin = 0; zbin < bins_Z; zbin++) {
                        double sp = sigma_total[0][zbin];
                        double sm = sigma_total[1][zbin];
                        if (sp == 0 || sm == 0) continue;
                        nz_filled++;

                        for (int c = 0; c < 2; c++) {
                            sum_dsig[c] += dsigma_sq_total[c][zbin];
                            for (int t = 0; t < 5; t++)
                                sum_term[c][t] += dsig_term[c][t][zbin];
                            for (int kk = 0; kk < 4; kk++)
                                sum_cat[c][kk] += dsig_cat[c][kk][zbin];
                        }
                        double cov = cov_total_z[zbin];
                        for (int t = 0; t < 4; t++)
                            sum_ct[t] += cov_term[t][zbin];

                        // Three master-formula contributions
                        double dsp_sq = dsigma_sq_total[0][zbin];
                        double dsm_sq = dsigma_sq_total[1][zbin];
                        sum_dRoR_sp  += dsp_sq / (sp * sp);
                        sum_dRoR_sm  += dsm_sq / (sm * sm);
                        sum_dRoR_cov += -2.0 * cov / (sp * sm);
                        sum_dRoR_sq  += dsp_sq / (sp * sp)
                                      + dsm_sq / (sm * sm)
                                      - 2.0 * cov / (sp * sm);
                    }

                    if (nz_filled > 0) {
                        const char* catname[4] = {
                            "counting stats  ",
                            "PID correction  ",
                            "VM weights      ",
                            "C_mc correction "
                        };
                        const char* tname[5] = {
                            "C_mc same-chg   ", "counting+PID-corr",
                            "VM same-chg     ", "C_mc cross-chg  ",
                            "VM cross-chg    "
                        };
                        const char* cname[4] = {
                            "Cov: -B+  * W+  * dC+^2",
                            "Cov: -B-  * W-  * dC-^2",
                            "Cov: +C+^2 * Var(W+)   ",
                            "Cov: +C-^2 * Var(W-)   "
                        };
                        cout << "\n──── Uncertainty breakdown: Q2 bin " << i
                             << ", xB bin " << j << ", var " << var
                             << " (" << nz_filled << " z-bins) ────\n";

                        // Four-source breakdown (user-requested)
                        for (int c = 0; c < 2; c++) {
                            cout << "  Var(sigma" << (c==0 ? "+" : "-")
                                 << ") — 4-source breakdown:\n";
                            if (sum_dsig[c] > 0)
                                for (int kk = 0; kk < 4; kk++)
                                    cout << "    " << catname[kk] << ":  "
                                         << Form("%6.2f%%", 100. * sum_cat[c][kk] / sum_dsig[c])
                                         << "\n";
                        }

                        // Full 5-term breakdown (for deeper debugging)
                        for (int c = 0; c < 2; c++) {
                            cout << "  Var(sigma" << (c==0 ? "+" : "-")
                                 << ") — 5-term detail:\n";
                            if (sum_dsig[c] > 0)
                                for (int t = 0; t < 5; t++)
                                    cout << "    " << tname[t] << ":  "
                                         << Form("%6.2f%%", 100. * sum_term[c][t] / sum_dsig[c])
                                         << "\n";
                        }

                        cout << "  Cov(sigma+, sigma-) breakdown:\n";
                        double sum_ct_abs = fabs(sum_ct[0]) + fabs(sum_ct[1])
                                          + fabs(sum_ct[2]) + fabs(sum_ct[3]);
                        if (sum_ct_abs > 0)
                            for (int t = 0; t < 4; t++)
                                cout << "    " << cname[t] << ":  "
                                     << Form("%+7.2f%%", 100. * sum_ct[t] / sum_ct_abs)
                                     << " (of |Cov| total)\n";

                        cout << "  (dR/R)^2 breakdown (summed over z):\n";
                        if (sum_dRoR_sq > 0) {
                            cout << "    (dsig+/sig+)^2 :  "
                                 << Form("%6.2f%%", 100. * sum_dRoR_sp / sum_dRoR_sq) << "\n";
                            cout << "    (dsig-/sig-)^2 :  "
                                 << Form("%6.2f%%", 100. * sum_dRoR_sm / sum_dRoR_sq) << "\n";
                            cout << "    -2*Cov/(s+*s-) :  "
                                 << Form("%+6.2f%%", 100. * sum_dRoR_cov / sum_dRoR_sq) << "\n";
                        }
                        cout << "\n";
                    }
                }

                // ── Write per-bin histograms ──
                outFile->cd();
                hZ      [var][i-1][j-1][0]->Write();
                hSigma  [var][i-1][j-1][0]->Write();
                hSigma  [var][i-1][j-1][1]->Write();
                hCovHist[var][i-1][j-1]   ->Write();
                hFF     [var][i-1][j-1]   ->Write();
            } // var
        } // j (xB)
    } // i (Q2)

    for (int E = 0; E < N_E; E++) delete corrByE[E];

    cout << "Closing Out File\n";
    delete outFile;
    cout << "Done\n";
    return 0;
}

// =============================================================================
// Variance-ingredient helper
// =============================================================================

// Var_PID_C^h = sum_p [ NMcPi^h(p)^2 * sigma_{C_pi2k^h(p)}^2
//                     + NMcK^h(p)^2  * sigma_{C_k2pi^h(p)}^2 ]
//
// Per-event generalisation of the reference's
// (C_mc^h)^2 * sum_p [ N_pi^h(p)^2 * sigma_{C_pi2k(p)}^2
//                    + N_K^h(p)^2  * sigma_{C_k2pi(p)}^2 ]
// where N_*(p) is replaced by NMc*(p) = sum of C_mc,i^h over events in
// the (p, charge) cell — so the per-event C_mc weighting is captured
// without ever evaluating C_mc at the bin centre.
static double computeVarPIDCorr_fullUnc(
    const std::vector<double>& NMcPi_p,
    const std::vector<double>& NMcK_p,
    const std::vector<double>& sigPi2k_p,
    const std::vector<double>& sigK2pi_p)
{
    double v = 0.0;
    int n = (int)NMcPi_p.size();
    for (int p = 0; p < n; p++){
        v += NMcPi_p[p] * NMcPi_p[p] * sigPi2k_p[p] * sigPi2k_p[p];
        v += NMcK_p [p] * NMcK_p [p] * sigK2pi_p[p] * sigK2pi_p[p];
    }
    return v;
}

// =============================================================================
// FF + uncertainty
// =============================================================================
//
//   R          = sigma+ / sigma-
//   Var(R)     = R^2 * [ Var(sigma+)/sigma+^2 + Var(sigma-)/sigma-^2
//                       - 2 Cov(sigma+, sigma-) / (sigma+ sigma-) ]
//   FF         = (4 - R) / (4R - 1)
//   dFF/dR     = -17 / (4R - 1)^2
//   Var(FF)    = ( 17 / (4R - 1)^2 )^2 * Var(R)
//   sigma_FF   = sqrt(Var(FF))
//
// Guards:
//   * sigma- == 0 or non-finite       -> valid=false, FF=NaN
//   * |4R - 1| < 0.05                 -> valid=true,  nearPole=true (warn)
//   * (dR/R)^2 negative due to roundoff -> clamp to 0 with warning
static FFResult_fullUnc computeFFAndUncertainty_fullUnc(
    double sp, double sm,
    double Var_sp, double Var_sm,
    double Cov_pm)
{
    FFResult_fullUnc out;
    out.R = std::nan(""); out.sigma_R = std::nan("");
    out.FF = std::nan(""); out.sigma_FF = std::nan("");
    out.valid = false; out.nearPole = false;

    const double kEps = 1e-12;
    if (!isfinite(sp) || !isfinite(sm) || fabs(sp) < kEps || fabs(sm) < kEps) return out;

    double R = sp / sm;
    double dRoR2 = Var_sp / (sp * sp) + Var_sm / (sm * sm)
                  - 2.0 * Cov_pm / (sp * sm);
    if (dRoR2 < 0.0) {
        cerr << "WARNING: (dR/R)^2 < 0 (numerical roundoff); clamping to 0.\n";
        dRoR2 = 0.0;
    }
    double Var_R = R * R * dRoR2;
    out.R = R;
    out.sigma_R = sqrt(fmax(Var_R, 0.0));

    double denom = 4.0 * R - 1.0;
    if (!isfinite(denom) || fabs(denom) < kEps) return out;
    if (fabs(denom) < 0.05) out.nearPole = true;

    out.FF = (4.0 - R) / denom;
    if (!isfinite(out.FF)) return out;

    double dFFdR = -17.0 / (denom * denom);
    out.sigma_FF = fabs(dFFdR) * out.sigma_R;
    out.valid = true;
    return out;
}

// =============================================================================
// Renamed legacy helpers (mirror originals; suffixed for symbol coexistence).
// =============================================================================

void zeroSuppress_fullUnc(TH1F * h){
    int nBins = h->GetNbinsX();
    for( int i = 1; i <= nBins; i++ ){
        if( h->GetBinContent(i) < 0 ){
            h->SetBinContent(i, 0);
            h->SetBinError(i, 0);
        }
    }
}

double getVarVal_fullUnc(TString var, electron e, pion pi){
    if( var == "null" ) return 0;
    if( var == "p_e" )      return e.get3Momentum().Mag();
    if( var == "theta_e" )  return e.get3Momentum().Theta()*rad_to_deg;
    if( var == "phi_e" )    return e.get3Momentum().Phi()*rad_to_deg;
    if( var == "W2" )       return e.getW2();
    if( var == "Q2" )       return e.getQ2();
    if( var == "xB" )       return e.getXb();
    if( var == "y" )        return e.getY();
    if( var == "sector" )   return e.getDC_sector();
    if( var == "sector_e" ){
        double phi_deg = e.get3Momentum().Phi()*rad_to_deg + 20.;
        if (phi_deg < 0.) phi_deg += 360.;
        return (int)(phi_deg / 60.) + 1;
    }
    if( var == "p_pi" )     return pi.get3Momentum().Mag();
    if( var == "theta_pi" ) return pi.get3Momentum().Theta()*rad_to_deg;
    if( var == "phi_pi" )   return pi.get3Momentum().Phi()*rad_to_deg;
    if( var == "phi_q" ){
        return pi.getPi_q().Phi()*rad_to_deg + 360*( (int)(pi.getPi_q().Phi() < 0) );
    }
    if( var == "Z" || var == "z" )    return pi.getZ();
    if( var == "Mx" || var == "M_x" ) return pi.getMx();
    if( var == "pT" || var == "Pt" )  return pi.getPi_q().Pt();
    if( var == "sector_pi" )          return pi.getDC_sector();
    if( var == "eta" )                return pi.getEta();
    return 0;
}

bool updateCorrectionsForBeam_fullUnc(
    double eBeam,
    double& beam_energy,
    TString& corrName,
    correctionTools& corrector)
{
    if (eBeam == beam_energy) return false;
    getNewRhoName_fullUnc(eBeam, corrName);
    corrector.setWeightName( corrName );
    if( n4d_corr_bins_fU > 0 ){
        corrector.setWeightName4D( corrName );
        corrector.setN4dBins(n4d_corr_bins_fU);
    }
    TString fitName = corrName;
    getNewRhoName_fullUnc(eBeam, fitName);
    // Only reload if the fit file name actually changed — avoids re-parsing
    // TF3 formulas (expensive) when corrections don't depend on beam energy.
    if( fitName != corrector.getWeightFitName() ){
        corrector.setWeightFitName( fitName );
        corrector.loadContinuousCorrections();
    }
    beam_energy = eBeam;
    return true;
}

bool updateCorrectionsForBeam_fullUnc(
    int E,
    TString& corrName,
    correctionTools& corrector)
{
    double eBeam = 10.2 + 0.2 * E;
    getNewRhoName_fullUnc(eBeam, corrName);
    corrector.setWeightName( corrName );
    if( n4d_corr_bins_fU > 0 ){
        corrector.setWeightName4D( corrName );
        corrector.setN4dBins(n4d_corr_bins_fU);
    }
    TString fitName = corrName;
    getNewRhoName_fullUnc(eBeam, fitName);
    corrector.setWeightFitName( fitName );
    corrector.loadContinuousCorrections();
    return true;
}

void getNewRhoName_fullUnc( double E, TString& rho_norm_name) {
    std::string fname = rho_norm_name.Data();
    std::ostringstream ss;
    ss << E;
    std::string newEnergyStr = ss.str();
    std::regex energyPattern(R"(10(\.\d+)?)");
    fname = std::regex_replace(fname, energyPattern, newEnergyStr,
                               std::regex_constants::format_first_only);
    rho_norm_name = TString(fname);
}
