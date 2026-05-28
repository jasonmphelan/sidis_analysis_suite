#include <cstdlib>
#include <iostream>
#include <vector>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TError.h"
#include "genElectron.h"
#include "genPion.h"
#include "analyzer.h"
#include "constants.h"

double getVarVal(TString var, genElectron e, genPion pi){
    if( var == "null" )      return 0;
    if( var == "Q2" )        return e.getQ2();
    if( var == "xB" )        return e.getXb();
    if( var == "W2" )        return e.getW2();
    if( var == "y" )         return e.getY();
    if( var == "p_e" )       return e.get3Momentum().Mag();
    if( var == "theta_e" )   return e.get3Momentum().Theta() * rad_to_deg;
    if( var == "phi_e" )     return e.get3Momentum().Phi()  * rad_to_deg;
    if( var == "Z"  || var == "z"   ) return pi.getZ();
    if( var == "Mx" || var == "M_x" ) return pi.getMx();
    if( var == "xF" )        return pi.getXf();
    if( var == "eta" )       return pi.getEta();
    if( var == "pT" || var == "Pt" ) return pi.getPi_q().Pt();
    if( var == "p_pi" )      return pi.get3Momentum().Mag();
    if( var == "theta_pi" )  return pi.get3Momentum().Theta() * rad_to_deg;
    if( var == "phi_pi" )    return pi.get3Momentum().Phi()   * rad_to_deg;
    return 0;
}

using std::cout;
using std::cerr;

// Usage: makeGenFFRatio [outFile] [target] [matchType] [nGen] [gen files...]
//          (opt: matchingFile) (opt: var_name bins_var var_min var_max) (opt: overflowBinning)
//   target:          0 = RGB/deuterium, 1 = RGA/proton
//   matchType:       0 = none, 1/2 = 2D matching, 3 = 2D+map (both charges), 4 = 2D+map (per charge)
//   var_name:        extra binning variable, e.g. pT, Mx, xF, eta, p_pi, theta_pi, Q2, y, W2
//   overflowBinning: bitmap — 1=Q2, 2=xB, 8=var; set bit to clamp out-of-range events into
//                    last bin and fill all bins up to the native bin (last bin = full integral)

int main(int argc, char** argv){

    if( argc < 5 ){
        cerr << "Usage: makeGenFFRatio [outFile] [target] [matchType] [nGen] [gen files...]\n";
        cerr << "       (opt: matchingFile) (opt: var_name bins_var var_min var_max) (opt: overflowBinning)\n";
        cerr << "  target:          0 = RGB/deuterium, 1 = RGA/proton\n";
        cerr << "  matchType:       0 = none, 1/2 = 2D, 3 = 2D+map (both), 4 = 2D+map (per charge)\n";
        cerr << "  var_name:        pT, Mx, xF, eta, p_pi, theta_pi, Q2, y, W2, etc.\n";
        cerr << "  overflowBinning: bitmap 1=Q2, 2=xB, 8=var (last bin = full integral)\n";
        return -1;
    }

    TString outName   = argv[1];
    int     target    = atoi(argv[2]);
    int     matchType = atoi(argv[3]);
    int     nGen      = atoi(argv[4]);
    const int firstGen = 5;

    if( argc < firstGen + nGen ){
        cerr << "Expected " << nGen << " gen files but only " << (argc - firstGen) << " provided.\n";
        return -1;
    }

    const int firstOpt = firstGen + nGen;
    TString matchingFile  = "matchCut2D_10.root";
    TString var_name      = "null";
    int     bins_var      = 1;
    double  var_min       = 0;
    double  var_max       = 1;
    int     overflowBinning = 0;  // bitmap: 1=Q2, 2=xB, 8=var
    if( argc > firstOpt )     matchingFile    = argv[firstOpt];
    if( argc > firstOpt + 1 ) var_name        = argv[firstOpt + 1];
    if( argc > firstOpt + 2 ) bins_var        = atoi(argv[firstOpt + 2]);
    if( argc > firstOpt + 3 ) var_min         = atof(argv[firstOpt + 3]);
    if( argc > firstOpt + 4 ) var_max         = atof(argv[firstOpt + 4]);
    if( argc > firstOpt + 5 ) overflowBinning = atoi(argv[firstOpt + 5]);
    int nVar = bins_var;

    TChain* genChain = new TChain("ePi");
    cerr << "Gen files:\n";
    for( int i = 0; i < nGen; i++ ){
        genChain->Add( argv[firstGen + i] );
        cerr << "  " << argv[firstGen + i] << "\n";
    }
    cerr << "matchType=" << matchType << "  matchingFile=" << matchingFile << "\n";
    if( nVar > 1 )
        cerr << "Extra binning: " << var_name << " [" << var_min << ", " << var_max << "] (" << nVar << " bins)\n";
    if( overflowBinning )
        cerr << "Overflow binning bitmap=" << overflowBinning
             << " (" << ((overflowBinning & 1) ? "Q2 " : "")
             << ((overflowBinning & 2) ? "xB " : "")
             << ((overflowBinning & 8) ? "var" : "") << ")\n";

    // ── Analyzer ───────────────────────────────────────────────────────────────
    analyzer anal(0, -1);
    anal.setAnalyzerLevel(1);
    anal.setTarget(0);
    anal.loadCutValues(-1, 10.2);
    anal.loadMatchingFunctions(matchingFile);
    anal.loadMatchingFunctions3D();
    anal.loadAcceptanceMapContinuous( (TString)_DATA + "/acceptance_map/acceptanceMap_allE_final.root" );

    // ── Histograms ─────────────────────────────────────────────────────────────
    const int nXb = bins_xB;
    const int nQ2 = bins_Q2;
    const int nZ  = bins_Z;

    // Index helper: (ix, iq, iv) -> flat index
    auto idx = [&](int ix, int iq, int iv){ return iv*nXb*nQ2 + ix*nQ2 + iq; };
    int nHist = nXb * nQ2 * nVar;

    std::vector<TH1F*> hPlus (nHist);
    std::vector<TH1F*> hMinus(nHist);
    std::vector<TH1F*> hRatio(nHist);
    std::vector<TH1F*> hFF   (nHist);
    std::vector<TH1F*> hDuP  (nHist);
    std::vector<TH1F*> hDuM  (nHist);
    std::vector<TH1F*> hDdP  (nHist);
    std::vector<TH1F*> hDdM  (nHist);

    for( int iv = 0; iv < nVar; iv++ ){
        for( int ix = 0; ix < nXb; ix++ ){
            for( int iq = 0; iq < nQ2; iq++ ){
                TString tag = Form("_%d_%d", iq+1, ix+1);
                if( nVar > 1 ) tag += Form("_v%d", iv+1);
                int i = idx(ix,iq,iv);
                hPlus [i] = new TH1F("hGenPlus"  + tag, ";z;N_{#pi^{+}}",               nZ, 0.3, 1.0);
                hMinus[i] = new TH1F("hGenMinus" + tag, ";z;N_{#pi^{-}}",               nZ, 0.3, 1.0);
                hRatio[i] = new TH1F("hRatio_0"  + tag, ";z;N_{#pi^{+}}/N_{#pi^{-}}",  nZ, 0.3, 1.0);
                hFF   [i] = new TH1F("hFF_0"     + tag, ";z;(4-r)/(4r-1)",              nZ, 0.3, 1.0);
                hDuP  [i] = new TH1F("hDuP_0"   + tag, ";z;D_{u}^{+}",                 nZ, 0.3, 1.0);
                hDuM  [i] = new TH1F("hDuM_0"   + tag, ";z;D_{u}^{-}",                 nZ, 0.3, 1.0);
                hDdP  [i] = new TH1F("hDdP_0"   + tag, ";z;D_{d}^{+}",                 nZ, 0.3, 1.0);
                hDdM  [i] = new TH1F("hDdM_0"   + tag, ";z;D_{d}^{-}",                 nZ, 0.3, 1.0);
                hPlus [i]->Sumw2();
                hMinus[i]->Sumw2();
                hFF   [i]->Sumw2();
            }
        }
    }

    // 2D overview: (xB, z) sliced by Q2 bin (integrated over var)
    std::vector<TH2F*> hFF2D(nQ2);
    for( int iq = 0; iq < nQ2; iq++ ){
        hFF2D[iq] = new TH2F(Form("hFF2D_%d", iq),
                             Form(";x_{B};z;(4-r)/(4r-1) (Q^{2} bin %d)", iq),
                             nXb, xB_min, xB_max, nZ, 0.3, 1.0);
    }

    // ── Event loop ─────────────────────────────────────────────────────────────
    TTreeReader reader(genChain);
    TTreeReaderValue<genElectron>  e_gen (reader, "e_gen");
    TTreeReaderArray<genPion>      pi_gen(reader, "pi_gen");
    TTreeReaderValue<int>      tar_string(reader, "string");

    int event_total = genChain->GetEntries();
    cout << "Total gen events: " << event_total << "\n";

    gErrorIgnoreLevel = kError;
    while( reader.Next() ){
        int ev = reader.GetCurrentEntry();
        if( ev % 1000000 == 0 )
            cout << "Events: " << ev << " / " << event_total << "\n";

        if( !anal.applyElectronKinematicCuts( *e_gen ) ) continue;

        double Q2 = e_gen->getQ2();
        double xB = e_gen->getXb();

        int bQ2 = (int)( ((Q2 - Q2_min) / (Q2_max - Q2_min)) * nQ2 );
        int bxB = (int)( ((xB - xB_min) / (xB_max - xB_min)) * nXb );

        // Reject under-range; clamp or reject over-range per overflow bitmap
        if( bQ2 < 0 || (!(overflowBinning & 1) && bQ2 >= nQ2) ) continue;
        if( bxB < 0 || (!(overflowBinning & 2) && bxB >= nXb) ) continue;
        if( (overflowBinning & 1) && bQ2 >= nQ2 ) bQ2 = nQ2 - 1;
        if( (overflowBinning & 2) && bxB >= nXb ) bxB = nXb - 1;

        for( auto pi : pi_gen ){
            double charge = pi.getCharge();
            int chargeIdx = (int)( charge < 1 );  // 0=pi+, 1=pi-


            int bVar = 0;
            if( nVar > 1 ){
                bVar = (int)(((getVarVal(var_name, *e_gen, pi) - var_min) / (var_max - var_min)) * nVar);
                if( bVar < 0 || (!(overflowBinning & 8) && bVar >= nVar) ) continue;
                if( (overflowBinning & 8) && bVar >= nVar ) bVar = nVar - 1;
            }

            // Overflow fill ranges: each event fills all bins from lo up to its native bin.
            // With one bin (integration mode) this is always just bin 0.
            int Q2_lo = (overflowBinning & 1) ? 0 : bQ2;
            int xB_lo = (overflowBinning & 2) ? 0 : bxB;
            int V_lo  = (overflowBinning & 8) ? 0 : bVar;


            //if( !anal.applyPionKinematicCuts(pi) || pi.getMx() < 1.7 ) continue;
            double phi    = pi.get3Momentum().Phi();
            double theta  = pi.get3Momentum().Theta() * rad_to_deg;
            double p      = pi.get3Momentum().Mag();

            // ── Sector from azimuthal angle (matches gen loop in calcCorrections) ──
            double sector_i = -1;
            if( charge > 0 ){
                if     ( phi >  -0.8  && phi <   0.25 ) sector_i = 1;
                else if( phi >=  0.25 && phi <   1.3  ) sector_i = 2;
                else if( phi >=  1.3  && phi <=  2.35 ) sector_i = 3;
                else if( phi >   2.35 || phi <  -2.9  ) sector_i = 4;
                else if( phi >  -2.9  && phi <  -1.85 ) sector_i = 5;
                else                                     sector_i = 6;
            } else {
                if     ( phi >  -0.25 && phi <   0.8  ) sector_i = 1;
                else if( phi >=  0.8  && phi <   1.85 ) sector_i = 2;
                else if( phi >=  1.85 && phi <=  2.9  ) sector_i = 3;
                else if( phi >   2.9  || phi <  -2.4  ) sector_i = 4;
                else if( phi >  -2.4  && phi <  -1.25 ) sector_i = 5;
                else                                     sector_i = 6;
            }
            
            if( !anal.applyPionKinematicCuts(pi) ) continue;


            // ── Acceptance matching ────────────────────────────────────────────
            bool matching = true;
            if( matchType == 1 || matchType == 2 ){
                matching = anal.acceptance_match_2d( theta, p, sector_i );
            } else if( matchType == 3 ){
                matching = anal.acceptance_match_2d( theta, p, sector_i ) &&
                           anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 1 ) >= 0 &&
                           anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, 2 ) >= 0;
            } else if( matchType == 4 ){
                matching = anal.acceptance_match_2d( theta, p, sector_i ) &&
                           anal.applyAcceptanceMap( p, rad_to_deg*phi, theta, chargeIdx + 1 ) >= 0;
            }
            if( !matching ) continue;
            
            double z = pi.getZ();
            if(*tar_string == 1){
                for(int iq=Q2_lo; iq<=bQ2; iq++) for(int ix=xB_lo; ix<=bxB; ix++) for(int iv=V_lo; iv<=bVar; iv++){
                    if( chargeIdx == 0 ) hDdP[idx(ix,iq,iv)]->Fill(pi.getZ());
                    if( chargeIdx == 1 ) hDdM[idx(ix,iq,iv)]->Fill(pi.getZ());
                }
            }
            if(*tar_string == 2){
                for(int iq=Q2_lo; iq<=bQ2; iq++) for(int ix=xB_lo; ix<=bxB; ix++) for(int iv=V_lo; iv<=bVar; iv++){
                    if( chargeIdx == 0 ) hDuP[idx(ix,iq,iv)]->Fill(pi.getZ());
                    if( chargeIdx == 1 ) hDuM[idx(ix,iq,iv)]->Fill(pi.getZ());
                }
            }

            for(int iq=Q2_lo; iq<=bQ2; iq++) for(int ix=xB_lo; ix<=bxB; ix++) for(int iv=V_lo; iv<=bVar; iv++){
                if( chargeIdx == 0 ) hPlus [idx(ix,iq,iv)]->Fill(z);
                else                 hMinus[idx(ix,iq,iv)]->Fill(z);
            }
        }
    }

    // ── Compute ratios and FF ──────────────────────────────────────────────────
    for( int iv = 0; iv < nVar; iv++ ){
        for( int ix = 0; ix < nXb; ix++ ){
            for( int iq = 0; iq < nQ2; iq++ ){
                int i = idx(ix,iq,iv);
                hRatio[i]->Divide( hPlus[i], hMinus[i], 1, 1, "B" );

                for( int iz = 1; iz <= nZ; iz++ ){
                    double r     = hRatio[i]->GetBinContent(iz);
                    double r_err = hRatio[i]->GetBinError(iz);

                    if( r > 0.25 ){  // (4r-1) > 0 required for physical FF
                        double denom  = 4.*r - 1.;
                        double ff     = (4. - r) / denom;
                        // d(FF)/dr = -15 / (4r-1)^2
                        double ff_err = 15. / (denom * denom) * r_err;
                        hFF[i]->SetBinContent(iz, ff);
                        hFF[i]->SetBinError  (iz, ff_err);
                        if( iv == 0 ){  // 2D summary only for first var bin
                            hFF2D[iq]->SetBinContent(ix + 1, iz, ff);
                            hFF2D[iq]->SetBinError  (ix + 1, iz, ff_err);
                        }
                    }
                }
            }
        }
    }

    // ── Write output ───────────────────────────────────────────────────────────
    TFile* outFile = new TFile(outName, "RECREATE");
    for( int iv = 0; iv < nVar; iv++ ){
        for( int ix = 0; ix < nXb; ix++ ){
            for( int iq = 0; iq < nQ2; iq++ ){
                int i = idx(ix,iq,iv);
                hPlus [i]->Write();
                hMinus[i]->Write();
                hRatio[i]->Write();
                hFF   [i]->Write();

                hDdP[i]->Sumw2();
                hDdM[i]->Sumw2();
                hDuP[i]->Sumw2();
                hDuM[i]->Sumw2();


                hDdP[i]->Divide(hDdM[i]);
                hDuM[i]->Divide(hDuP[i]);

                hDdP[i]->Write();
                hDuP[i]->Write();
                hDdM[i]->Write();
                hDuM[i]->Write();
            }
        }
    }
    for( int iq = 0; iq < nQ2; iq++ )
        hFF2D[iq]->Write();

    outFile->Close();
    delete outFile;
    delete genChain;
    cout << "Done. Output: " << outName << "\n";
    return 0;
}
