
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMultiGraph.h"

using std::cout;
using std::cerr;

// Acceptance map parameters indexed as [particle][sector][param]
// params: 0=lower_a, 1=upper_a, 2=lower_b, 3=upper_b, 4=theta_min, 5=theta_max, 6=phi_mean
TF1 * mapParameters[3][6][7];

double mapFunc(int particle, int sec, int param, double p){
    if( mapParameters[particle][sec][param] == nullptr ) return 0;
    return mapParameters[particle][sec][param]->Eval(p);
}

// Returns phi_min and phi_max at a given (theta, p) for a given sector+particle.
// phi_min and phi_max are returned in the "raw" frame before wrapping.
bool getPhiBounds(int particle, int sec, double p, double theta,
                  double &phi_min, double &phi_max){

    double map_min = mapFunc(particle, sec, 4, p);
    double map_max = mapFunc(particle, sec, 5, p);
    if( theta <= map_min || theta >= map_max ) return false;

    double phi_avg = mapFunc(particle, sec, 6, p);
    double a_low   = mapFunc(particle, sec, 0, p);
    double a_up    = mapFunc(particle, sec, 1, p);
    double b_low   = mapFunc(particle, sec, 2, p);
    double b_up    = mapFunc(particle, sec, 3, p);

    double denom_up, denom_low;
    if( particle == 1 ){
        denom_up  = exp( (theta - map_min)/b_up  );
        denom_low = exp( (theta - map_min)/b_low );
    } else {
        denom_up  = (theta - map_min)/b_up;
        denom_low = (theta - map_min)/b_low;
    }

    phi_max = phi_avg + a_up  * (1. - 1./( denom_up  + 1. ));
    phi_min = phi_avg - a_low * (1. - 1./( denom_low + 1. ));
    return true;
}

int main( int argc, char** argv ){

    if( argc < 3 ){
        cerr << "Usage: ./plotAcceptanceRegions [acceptance_map.root] [output.root] [p1] [p2] ...\n";
        cerr << "  If no momenta are given, defaults: 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 GeV/c\n";
        return -1;
    }

    TString mapFile  = argv[1];
    TString outFile  = argv[2];

    // Collect momentum values from argv or use defaults
    std::vector<double> pVals;
    if( argc > 3 ){
        for( int i = 3; i < argc; i++ )
            pVals.push_back( atof(argv[i]) );
    } else {
        pVals = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    }

    // Load acceptance map
    TFile * fMap = TFile::Open(mapFile);
    if( !fMap || fMap->IsZombie() ){
        cerr << "Cannot open acceptance map file: " << mapFile << "\n";
        return -1;
    }

    TString parName[3] = {"e", "pip", "pim"};
    TString paramName[7] = {"lower_a","upper_a","lower_b","upper_b","theta","max","mean"};

    for( int par = 0; par < 3; par++ ){
        for( int sec = 0; sec < 6; sec++ ){
            for( int param = 0; param < 7; param++ ){
                TString key = Form("f_") + parName[par] + Form("_%i_param_", sec+1) + paramName[param];
                mapParameters[par][sec][param] = (TF1*) fMap->Get( key );
                if( !mapParameters[par][sec][param] ){
                    cerr << "Warning: could not find " << key << "\n";
                }
            }
        }
    }

    TFile * outf = new TFile( outFile, "RECREATE" );

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    // Sector colors (one per sector)
    int secColor[6] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1, kCyan+2};

    // Particle labels for display
    TString parLabel[3] = {"Electron", "#pi^{+}", "#pi^{-}"};

    // Number of theta points for the boundary curves
    const int nThetaPts = 500;

    // For each particle type, draw all sectors and momenta on one canvas per p
    for( int par = 0; par < 3; par++ ){

        for( double p : pVals ){

            TString cName = Form("c_%s_p%.2f", parName[par].Data(), p);
            TCanvas * c = new TCanvas( cName, cName, 1200, 800 );

            // Find global theta range for this particle at this p
            double theta_min_all =  999.;
            double theta_max_all = -999.;
            for( int sec = 0; sec < 6; sec++ ){
                double tmin = mapFunc(par, sec, 4, p);
                double tmax = mapFunc(par, sec, 5, p);
                if( tmin < theta_min_all ) theta_min_all = tmin;
                if( tmax > theta_max_all ) theta_max_all = tmax;
            }
            // Add a small margin
            theta_min_all = std::max(0., theta_min_all - 1.);
            theta_max_all = theta_max_all + 1.;

            // Phi range: -40 to 360 covers the full rotated range
            TH2F * hFrame = new TH2F( Form("hFrame_%s_p%.2f", parName[par].Data(), p),
                                      Form("%s acceptance regions, p = %.2f GeV/c;#phi (deg);#theta (deg)",
                                           parLabel[par].Data(), p),
                                      1, -40, 370, 1, theta_min_all, theta_max_all );
            hFrame->Draw();

            TLegend * leg = new TLegend(0.75, 0.65, 0.92, 0.92);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.03);

            for( int sec = 0; sec < 6; sec++ ){

                double tmin = mapFunc(par, sec, 4, p);
                double tmax = mapFunc(par, sec, 5, p);
                if( tmin >= tmax ) continue;

                std::vector<double> phi_lo_pts, phi_hi_pts, theta_pts_lo, theta_pts_hi;

                for( int ti = 0; ti < nThetaPts; ti++ ){
                    double theta = tmin + (tmax - tmin) * ti / (nThetaPts - 1.);
                    double phi_min, phi_max;
                    if( !getPhiBounds(par, sec, p, theta, phi_min, phi_max) ) continue;

                    // Apply same phi wrapping logic as applyAcceptanceMap
                    // (here we just draw in the raw frame; the map stores absolute phi)
                    phi_lo_pts.push_back( phi_min );
                    phi_hi_pts.push_back( phi_max );
                    theta_pts_lo.push_back( theta );
                    theta_pts_hi.push_back( theta );
                }

                if( phi_lo_pts.empty() ) continue;
                int n = phi_lo_pts.size();

                // Build a closed polygon: lower boundary forward, upper boundary reversed
                std::vector<double> poly_phi, poly_theta;
                for( int i = 0; i < n; i++ ){
                    poly_phi.push_back( phi_lo_pts[i] );
                    poly_theta.push_back( theta_pts_lo[i] );
                }
                for( int i = n-1; i >= 0; i-- ){
                    poly_phi.push_back( phi_hi_pts[i] );
                    poly_theta.push_back( theta_pts_hi[i] );
                }
                // Close polygon
                poly_phi.push_back( phi_lo_pts[0] );
                poly_theta.push_back( theta_pts_lo[0] );

                TGraph * gPoly = new TGraph( poly_phi.size(), poly_phi.data(), poly_theta.data() );
                gPoly->SetName( Form("gAcc_%s_sec%i_p%.2f", parName[par].Data(), sec+1, p) );
                gPoly->SetLineColor( secColor[sec] );
                gPoly->SetLineWidth(2);
                gPoly->SetFillColor( secColor[sec] );
                gPoly->SetFillStyle(3001);
                gPoly->Draw("F SAME");
                gPoly->Draw("L SAME");

                // Draw theta_min and theta_max lines only for the phi span
                double phi_span_min = *std::min_element(phi_lo_pts.begin(), phi_lo_pts.end());
                double phi_span_max = *std::max_element(phi_hi_pts.begin(), phi_hi_pts.end());
                TLine * lMin = new TLine( phi_span_min, tmin, phi_span_max, tmin );
                TLine * lMax = new TLine( phi_span_min, tmax, phi_span_max, tmax );
                lMin->SetLineColor( secColor[sec] );
                lMax->SetLineColor( secColor[sec] );
                lMin->SetLineWidth(2);
                lMax->SetLineWidth(2);
                lMin->Draw("SAME");
                lMax->Draw("SAME");

                if( par == 0 ){
                    leg->AddEntry( gPoly, Form("Sector %i", sec+1), "f" );
                } else {
                    leg->AddEntry( gPoly, Form("Sec %i", sec+1), "f" );
                }

                gPoly->Write();
            }

            leg->Draw();

            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.04);
            lat.DrawLatex(0.12, 0.92, Form("%s, p = %.2f GeV/c", parLabel[par].Data(), p));

            c->Write();
            c->SaveAs( Form("accRegions_%s_p%.2f.pdf", parName[par].Data(), p) );

            delete c;
        }

        // Also make a single canvas per particle showing all p values overlaid for one representative sector
        // (sector 0, i.e. sector 1)
        int refSec = 0;
        TCanvas * cOverlay = new TCanvas( Form("c_%s_overlay_sec1", parName[par].Data()),
                                          Form("%s sector 1 overlay all p", parName[par].Data()),
                                          900, 700 );

        TH2F * hOvFrame = new TH2F( Form("hOvFrame_%s", parName[par].Data()),
                                     Form("%s acceptance region, sector 1;#phi (deg);#theta (deg)",
                                          parLabel[par].Data()),
                                     1, -40, 370, 1, 0, 45 );
        hOvFrame->Draw();

        TLegend * legOv = new TLegend(0.75, 0.55, 0.92, 0.92);
        legOv->SetBorderSize(0);
        legOv->SetTextSize(0.03);

        int pIdx = 0;
        int pColors[] = {kRed, kOrange+1, kYellow+2, kGreen+2, kCyan+2, kBlue+1, kViolet+1, kMagenta+1};
        for( double p : pVals ){
            double tmin = mapFunc(par, refSec, 4, p);
            double tmax = mapFunc(par, refSec, 5, p);
            if( tmin >= tmax ){ pIdx++; continue; }

            std::vector<double> phi_lo, phi_hi, th_lo, th_hi;
            for( int ti = 0; ti < nThetaPts; ti++ ){
                double theta = tmin + (tmax - tmin)*ti/(nThetaPts-1.);
                double phi_min, phi_max;
                if( !getPhiBounds(par, refSec, p, theta, phi_min, phi_max) ) continue;
                phi_lo.push_back(phi_min);
                phi_hi.push_back(phi_max);
                th_lo.push_back(theta);
                th_hi.push_back(theta);
            }
            if( phi_lo.empty() ){ pIdx++; continue; }
            int n = phi_lo.size();

            std::vector<double> poly_phi, poly_theta;
            for(int i=0;i<n;i++){ poly_phi.push_back(phi_lo[i]); poly_theta.push_back(th_lo[i]); }
            for(int i=n-1;i>=0;i--){ poly_phi.push_back(phi_hi[i]); poly_theta.push_back(th_hi[i]); }
            poly_phi.push_back(phi_lo[0]); poly_theta.push_back(th_lo[0]);

            TGraph * g = new TGraph(poly_phi.size(), poly_phi.data(), poly_theta.data());
            g->SetLineColor( pColors[pIdx % 8] );
            g->SetLineWidth(2);
            g->Draw("L SAME");

            legOv->AddEntry(g, Form("p = %.2f GeV/c", p), "l");
            g->Write( Form("gOvAcc_%s_sec1_p%.2f", parName[par].Data(), p) );
            pIdx++;
        }
        legOv->Draw();
        cOverlay->Write();
        cOverlay->SaveAs( Form("accRegions_%s_overlay_sec1.pdf", parName[par].Data()) );
        delete cOverlay;
    }

    outf->Close();
    fMap->Close();
    cout << "Done. Output written to " << outFile << "\n";
    return 0;
}
