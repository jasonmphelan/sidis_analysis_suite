// fixRatioErrors.C
// Recomputes hRatio_p and hRatio_theta_p histograms from kRICH output files
// using binomial error propagation (hP_eb is a strict subset of hP_true).
//
// Usage:
//   root -l -b -q 'fixRatioErrors.C("kRICH_final.root")'
//   root -l -b -q 'fixRatioErrors.C("kRICH_final.root", "kRICH_fixed.root")'
//
// If no output name is given, "_fixed" is appended to the input name.

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>

void calcFactor( TH1F * piFactor, TH1F * piRatio, TH1F * kRatio );
void calcKaonFactor( TH1F * kFactor, TH1F * piRatio, TH1F * kRatio );

void fixRatioErrors( TString inName, TString outName = "" ){

    if( outName == "" ){
        outName = inName;
        outName.ReplaceAll(".root", "_fixed.root");
    }

    TFile * inFile = TFile::Open( inName, "READ" );
    if( !inFile || inFile->IsZombie() ){
        std::cerr << "Could not open input file: " << inName << "\n";
        return;
    }

    TFile * outFile = new TFile( outName, "RECREATE" );

    // species names must match what richKaonCorrections.cpp wrote
    // hP_eb / hP_true use the "species" array: pip, pim, kp, km, prot, aprot
    // ratio histograms are named with "ratio_species": pip, pim, kp, km, prtp, prtm
    const int N = 6;
    TString species[N]       = {"pip", "pim", "kp", "km", "prot", "aprot"};
    TString ratio_species[N] = {"pip", "pim", "kp", "km", "prtp", "prtm"};

    TString titles_1d[N] = {
        "#pi^{+} ID efficiency;p (GeV/c);#epsilon",
        "#pi^{-} ID efficiency;p (GeV/c);#epsilon",
        "K^{+} ID efficiency;p (GeV/c);#epsilon",
        "K^{-} ID efficiency;p (GeV/c);#epsilon",
        "p ID efficiency;p (GeV/c);#epsilon",
        "#bar{p} ID efficiency;p (GeV/c);#epsilon"
    };

    TH1F * hRatio_pip = nullptr;
    TH1F * hRatio_pim = nullptr;
    TH1F * hRatio_kp  = nullptr;
    TH1F * hRatio_km  = nullptr;

    for( int i = 0; i < N; i++ ){

        // --- 1D ratio ---
        TH1F * hNum = (TH1F*) inFile->Get( "hP_eb_"   + species[i] );
        TH1F * hDen = (TH1F*) inFile->Get( "hP_true_" + species[i] );

        if( !hNum || !hDen ){
            std::cerr << "Missing histograms for species " << species[i] << ", skipping.\n";
            continue;
        }

        // Clone into output file
        outFile->cd();
        TH1F * hRatio = (TH1F*) hNum->Clone( "hRatio_p_" + ratio_species[i] );
        hRatio->Reset();

        // "B" option: binomial errors  sigma = sqrt(eps*(1-eps)/N_true)
        hRatio->Divide( hNum, hDen, 1.0, 1.0, "B" );
        hRatio->SetTitle( titles_1d[i] );
        hRatio->Write();

        if( species[i] == "pip" ) hRatio_pip = hRatio;
        if( species[i] == "pim" ) hRatio_pim = hRatio;
        if( species[i] == "kp"  ) hRatio_kp  = hRatio;
        if( species[i] == "km"  ) hRatio_km  = hRatio;

        

        // --- 2D ratio (theta vs p) ---
        TH2F * hNum2 = (TH2F*) inFile->Get( "hTheta_p_eb_"   + species[i] );
        TH2F * hDen2 = (TH2F*) inFile->Get( "hTheta_p_true_" + species[i] );

        if( hNum2 && hDen2 ){
            TH2F * hRatio2 = (TH2F*) hNum2->Clone( "hRatio_theta_p_" + ratio_species[i] );
            hRatio2->Reset();
            hRatio2->Divide( hNum2, hDen2, 1.0, 1.0, "B" );
            hRatio2->SetTitle(
                "#theta vs p efficiency (" + ratio_species[i] + ");p (GeV/c);#theta (deg)" );
            hRatio2->Write();
        } else {
            std::cerr << "Missing 2D histograms for species " << species[i] << ", skipping 2D.\n";
        }

        // Also copy the raw eb/true histograms so the output is self-contained
        TH1F * hNum_copy = (TH1F*) hNum->Clone();
        TH1F * hDen_copy = (TH1F*) hDen->Clone();
        hNum_copy->Write();
        hDen_copy->Write();
    }

    // ── Correction factors ───────────────────────────────────────────────
    if( hRatio_pip && hRatio_kp ){
        TH1F * hFactor_pip = (TH1F*) hRatio_pip->Clone("hFactor_pip");
        hFactor_pip->Reset();
        hFactor_pip->SetTitle("#pi^{+} correction factor;p (GeV/c);correction");
        calcFactor( hFactor_pip, hRatio_pip, hRatio_kp );
        hFactor_pip->Write();
    } else {
        std::cerr << "Warning: missing pip or kp ratio — skipping hFactor_pip.\n";
    }

    if( hRatio_pim && hRatio_km ){
        TH1F * hFactor_pim = (TH1F*) hRatio_pim->Clone("hFactor_pim");
        hFactor_pim->Reset();
        hFactor_pim->SetTitle("#pi^{-} correction factor;p (GeV/c);correction");
        calcFactor( hFactor_pim, hRatio_pim, hRatio_km );
        hFactor_pim->Write();
    } else {
        std::cerr << "Warning: missing pim or km ratio — skipping hFactor_pim.\n";
    }

    if( hRatio_pip && hRatio_kp ){
        TH1F * hFactor_kp = (TH1F*) hRatio_kp->Clone("hFactor_kp");
        hFactor_kp->Reset();
        hFactor_kp->SetTitle("K^{+} correction factor;p (GeV/c);correction");
        calcKaonFactor( hFactor_kp, hRatio_pip, hRatio_kp );
        hFactor_kp->Write();
    } else {
        std::cerr << "Warning: missing pip or kp ratio — skipping hFactor_kp.\n";
    }

    if( hRatio_pim && hRatio_km ){
        TH1F * hFactor_km = (TH1F*) hRatio_km->Clone("hFactor_km");
        hFactor_km->Reset();
        hFactor_km->SetTitle("K^{-} correction factor;p (GeV/c);correction");
        calcKaonFactor( hFactor_km, hRatio_pim, hRatio_km );
        hFactor_km->Write();
    } else {
        std::cerr << "Warning: missing pim or km ratio — skipping hFactor_km.\n";
    }

    // Copy any other histograms from the input (nRICH, probRICH, etc.)
    TIter next( inFile->GetListOfKeys() );
    TKey * key;
    while( (key = (TKey*) next()) ){
        TString name = key->GetName();
        // Skip ones we already wrote
        if( name.BeginsWith("hRatio_") ) continue;
        if( name.BeginsWith("hP_eb_") )   continue;
        if( name.BeginsWith("hP_true_") ) continue;
        if( name.BeginsWith("hTheta_p_eb_") )   continue;
        if( name.BeginsWith("hTheta_p_true_") ) continue;

        TObject * obj = key->ReadObj();
        if( obj ){
            outFile->cd();
            obj->Write();
        }
    }

    outFile->Close();
    inFile->Close();

    std::cout << "Done. Output written to: " << outName << "\n";
}

void calcKaonFactor(TH1F * kFactor, TH1F * piRatio, TH1F * kRatio){

    int nBins = kFactor->GetNbinsX();
    for( int bin = 1; bin <= nBins; bin++ ){
        double cPi  = piRatio->GetBinContent(bin);
        double cK   = kRatio->GetBinContent(bin);
        double ePi  = piRatio->GetBinError(bin);
        double eK   = kRatio->GetBinError(bin);

        double denom  = cPi + cK - 1.;
        double corrK  = (1. - cK) / denom;
        kFactor->SetBinContent( bin, corrK );

        // df/d(cK)  = -cPi / denom^2
        // df/d(cPi) = -(1-cK) / denom^2
        if( std::abs(denom) > 0 ){
            double dfdK  = -cPi       / (denom * denom);
            double dfdPi = -(1. - cK) / (denom * denom);
            double err   = std::sqrt( dfdK * dfdK * eK * eK + dfdPi * dfdPi * ePi * ePi );
            kFactor->SetBinError( bin, err );
        } else {
            kFactor->SetBinError( bin, 0. );
        }
    }
}

void calcFactor(TH1F * piFactor, TH1F * piRatio, TH1F * kRatio){

    int nBins = piFactor->GetNbinsX();
    for ( int bin = 1 ; bin <= nBins ; bin++ ){
        double cPi  = piRatio->GetBinContent(bin);
        double cK   = kRatio->GetBinContent(bin);
        double ePi  = piRatio->GetBinError(bin);
        double eK   = kRatio->GetBinError(bin);

        double denom  = cK + cPi - 1.;
        double corrPi = cK / denom;
        piFactor->SetBinContent( bin, corrPi );

        // Error propagation: d(corrPi)/d(cK) = (cPi-1)/denom^2
        //                    d(corrPi)/d(cPi) = -cK/denom^2
        if( std::abs(denom) > 0 ){
            double dfdK  = (cPi - 1.) / (denom * denom);
            double dfdPi = -cK        / (denom * denom);
            double err   = std::sqrt( dfdK * dfdK * eK * eK + dfdPi * dfdPi * ePi * ePi );
            piFactor->SetBinError( bin, err );
        } else {
            piFactor->SetBinError( bin, 0. );
        }
    }

}
