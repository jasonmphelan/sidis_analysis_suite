// compare_Rd.C
// After all threshold scans are complete, read the per-threshold R_d files
// and print a side-by-side comparison table.
// Also compare to HERMES reference values.
//
// Usage: root -l -b -q 'compare_Rd.C'

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

void compare_Rd() {
    const int nThresholds = 5;
    int thresholds[nThresholds] = {1, 5, 10, 15, 20};

    const int nZ = 14;
    double zCenter[nZ];
    double rd[nThresholds][nZ];
    double rdErr[nThresholds][nZ];

    // Initialize
    for (int t = 0; t < nThresholds; t++)
        for (int z = 0; z < nZ; z++) {
            rd[t][z] = 0;
            rdErr[t][z] = 0;
        }

    // Read each file
    TString scanDir = gSystem->GetFromPipe("cd \"$(dirname \"$0\")\" 2>/dev/null || true; pwd");
    // Fallback: try to find files relative to current directory or standard location
    std::vector<TString> searchPaths;
    searchPaths.push_back("theta_scan/");
    searchPaths.push_back("../theta_scan/");
    searchPaths.push_back("");

    for (int t = 0; t < nThresholds; t++) {
        TString fname;
        bool found = false;
        for (auto& prefix : searchPaths) {
            fname = Form("%sRd_theta_%d.txt", prefix.Data(), thresholds[t]);
            std::ifstream test(fname.Data());
            if (test.good()) { found = true; test.close(); break; }
        }
        if (!found) {
            std::cerr << "WARNING: Could not find Rd_theta_" << thresholds[t] << ".txt" << std::endl;
            continue;
        }

        std::ifstream in(fname.Data());
        std::string line;
        std::getline(in, line); // skip header
        int z = 0;
        while (std::getline(in, line) && z < nZ) {
            std::istringstream iss(line);
            double zc, r, re;
            if (iss >> zc >> r >> re) {
                zCenter[z] = zc;
                rd[t][z] = r;
                rdErr[t][z] = re;
                z++;
            }
        }
        in.close();
    }

    // HERMES reference values (approximate)
    const int nRef = 5;
    double zRef[nRef]  = {0.32, 0.38, 0.44, 0.52, 0.62};
    double rdRef[nRef] = {1.68, 1.78, 1.90, 2.05, 2.22};

    // Print header
    std::cout << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << "  R_d = sigma(pi+)/sigma(pi-) vs z for each matching threshold" << std::endl;
    std::cout << "================================================================================" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(8) << "z";
    for (int t = 0; t < nThresholds; t++)
        std::cout << std::setw(12) << Form("%d%%", thresholds[t]);
    std::cout << std::setw(12) << "HERMES" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;

    for (int z = 0; z < nZ; z++) {
        std::cout << std::setw(8) << zCenter[z];
        for (int t = 0; t < nThresholds; t++) {
            if (rd[t][z] > 0)
                std::cout << std::setw(12) << rd[t][z];
            else
                std::cout << std::setw(12) << "---";
        }

        // Find closest HERMES reference
        double minDist = 999;
        int closest = -1;
        for (int r = 0; r < nRef; r++) {
            double d = fabs(zCenter[z] - zRef[r]);
            if (d < minDist && d < 0.04) { minDist = d; closest = r; }
        }
        if (closest >= 0)
            std::cout << std::setw(12) << rdRef[closest];
        else
            std::cout << std::setw(12) << "";
        std::cout << std::endl;
    }

    // Print mean R_d for z = 0.3-0.7
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    std::cout << std::setw(8) << "Mean";
    for (int t = 0; t < nThresholds; t++) {
        double sum = 0; int n = 0;
        for (int z = 0; z < nZ; z++) {
            if (zCenter[z] <= 0.7 && rd[t][z] > 0) { sum += rd[t][z]; n++; }
        }
        if (n > 0) std::cout << std::setw(12) << sum/n;
        else       std::cout << std::setw(12) << "---";
    }
    std::cout << std::endl;

    // Find which threshold gives closest match to HERMES
    std::cout << std::endl;
    std::cout << "=== Closest match to HERMES reference ===" << std::endl;
    for (int t = 0; t < nThresholds; t++) {
        double chi2 = 0; int nPts = 0;
        for (int r = 0; r < nRef; r++) {
            // Find z bin closest to reference
            int bestZ = -1; double bestDist = 999;
            for (int z = 0; z < nZ; z++) {
                double d = fabs(zCenter[z] - zRef[r]);
                if (d < bestDist) { bestDist = d; bestZ = z; }
            }
            if (bestZ >= 0 && rd[t][bestZ] > 0 && rdErr[t][bestZ] > 0) {
                double pull = (rd[t][bestZ] - rdRef[r]) / rdErr[t][bestZ];
                chi2 += pull * pull;
                nPts++;
            }
        }
        std::cout << Form("  %2d%%: chi2/ndf = %.2f / %d = %.2f",
                          thresholds[t], chi2, nPts, nPts > 0 ? chi2/nPts : 0) << std::endl;
    }
    std::cout << std::endl;
}
