// extract_Rd.C
// Extract R_d = sigma(pi+)/sigma(pi-) vs z from makeRatioBinned3D output.
// The output file contains hRatio_0_i_j histograms (var=0, xB bin i, Q2 bin j)
// where the pi+/pi- division has already been done.
// This macro averages R_d over all xB and Q2 bins (weighted by bin content != 0)
// and writes a text file with z_center and R_d columns.
//
// Usage: root -l -b -q 'extract_Rd.C("ratio_file.root", "output.txt")'

#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <cmath>

void extract_Rd(TString inFileName, TString outFileName) {
    TFile *f = TFile::Open(inFileName);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << inFileName << std::endl;
        return;
    }

    const int nXB = 14;
    const int nQ2 = 12;
    const int nZ  = 14;
    const double zMin = 0.3;
    const double zMax = 1.0;
    const double dz = (zMax - zMin) / nZ;

    // Accumulate weighted average of R_d across xB, Q2 bins
    double sumRd[nZ]    = {0};
    double sumWeight[nZ] = {0};
    int nHists = 0;
    int nWarnings = 0;

    for (int i = 1; i <= nXB; i++) {
        for (int j = 1; j <= nQ2; j++) {
            TString hname = Form("hRatio_0_%d_%d", i, j);
            TH1F *h = (TH1F*)f->Get(hname);
            if (!h) continue;

            // Check if histogram has any content
            bool hasContent = false;
            for (int k = 1; k <= nZ; k++) {
                if (h->GetBinContent(k) != 0) { hasContent = true; break; }
            }
            if (!hasContent) continue;
            nHists++;

            for (int k = 1; k <= nZ; k++) {
                double val = h->GetBinContent(k);
                double err = h->GetBinError(k);

                if (val == 0 || !std::isfinite(val)) continue;

                // Flag suspicious values
                if (val > 10 || val < 0) {
                    nWarnings++;
                    if (nWarnings <= 5)
                        std::cerr << "WARNING: R_d=" << val << " in " << hname
                                  << " z-bin " << k << std::endl;
                    continue;  // skip pathological bins
                }

                // Weight by 1/err^2 if error available, else equal weight
                double w = (err > 0) ? 1.0 / (err * err) : 1.0;
                sumRd[k-1]    += val * w;
                sumWeight[k-1] += w;
            }
        }
    }

    std::cout << "Found " << nHists << " non-empty ratio histograms" << std::endl;
    if (nWarnings > 0)
        std::cerr << "Total warnings (R_d out of range): " << nWarnings << std::endl;

    // Write output
    std::ofstream out(outFileName.Data());
    out << "z_center  Rd  Rd_err" << std::endl;
    for (int k = 0; k < nZ; k++) {
        double z = zMin + (k + 0.5) * dz;
        double rd = 0, rd_err = 0;
        if (sumWeight[k] > 0) {
            rd = sumRd[k] / sumWeight[k];
            rd_err = 1.0 / sqrt(sumWeight[k]);
        }
        out << Form("%.4f  %.4f  %.4f", z, rd, rd_err) << std::endl;
        std::cout << Form("z=%.3f  R_d=%.4f +/- %.4f", z, rd, rd_err) << std::endl;
    }
    out.close();
    std::cout << "Written to " << outFileName << std::endl;

    f->Close();
}
