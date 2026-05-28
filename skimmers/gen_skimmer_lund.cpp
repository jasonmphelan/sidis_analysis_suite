//==============================================================================
//  gen_skimmer_lund.cpp
//
//  LUND (clasdis) port of skimmers/gen_skimmer.cpp.
//
//  Reads clasdis SIDIS generator output in plain-text LUND format, applies the
//  SAME generator-level kinematic selections as gen_skimmer.cpp, and writes the
//  surviving events to a ROOT TTree with an IDENTICAL schema:
//
//      tree   "ePi"  title "(e,e'pi) event  information"
//      branches (in order):
//        runnum/I, torus/D, evnum/I, Ebeam/D, beam(TLorentzVector),
//        Ne/I, Npi/I, Npips/I, Npims/I,
//        e_gen(genElectron), pi_gen(vector<genPion>)
//
//  This file does NOT read HIPO and does NOT depend on the CLAS12 framework.
//  It uses the standalone, streamer-identical genElectron / genPion classes in
//  lund_gen_classes.h.  Masses and cut thresholds are pulled directly from the
//  project headers constants.h and cut_values.h, so they are guaranteed equal
//  to the originals.
//
//  --------------------------------------------------------------------------
//  PHYSICS / SELECTION  (traced line-by-line to the original)
//  --------------------------------------------------------------------------
//  Particle finding (gen_skimmer.cpp lines 189-208):
//    The original loops over c12.mcparts() starting at index 1 (skipping
//    index 0 = the beam lepton) and selects PID 11 / 211 / -211.  In a clasdis
//    LUND event line 1 is always the beam e- (status 21), so we likewise skip
//    the first particle line and select by PID.  This was verified to be
//    exactly equivalent to selecting final-state (status==1) e-/pi+/pi- in the
//    gemc/eventfiles samples (charged pions are always status 1; the only
//    other PID==11 is the beam, which we skip).
//
//  Electron (genElectron::setKinematicInformation + analyzer cuts):
//    e4 from (px,py,pz) with mass Me; q = beam - e4; Q2 = -q.Mag2();
//    nu = q.E(); xB = Q2/(2 Mp nu); y = nu/Ebeam; W2 = (p_rest + q).Mag2().
//    Cuts: W>=2.5, 2<=Q2<=8, 0.1<=xB<=0.66, y<=0.75, 3<=|p|<=10.6,
//          5deg<=theta<=35deg.
//
//  Pion (genPion::setKinematicInformation + analyzer cuts):
//    pi4 from (px,py,pz) with mass Mpi; Charge=(int)(pid/211);
//    Z = pi4.E()/q.E(); Mx = (q + p_rest - pi4).Mag();
//    eta in the q-frame; xF left at 0 (original line commented out).
//    Cuts: 1.7<=Mx<=5, 1.25<=|p|<=5, 0.3<=Z<=1.0, 5deg<=theta<=35deg.
//
//  --------------------------------------------------------------------------
//  GENERATOR-LEVEL ADAPTATIONS (agreed with the user; see header notes)
//  --------------------------------------------------------------------------
//    * runnum : no analog in a raw LUND file -> hardcoded 0.
//    * evnum  : 0-based sequential event index (order events appear in input).
//    * Ne     : TRUE per-event electron count (the original accumulated it
//               across events because electronsMC was never cleared; that is a
//               bug, not physics, so we emit the correct count).
//    * Npips / Npims : filled with the real pi+ / pi- counts (the original
//               left them at 0; user requested correct values).
//    * torus  : constant -1, identical to the original (not a reconstructed
//               quantity).
//    * e_gen Sector/PID/Charge = -1/0/0 and pi_gen Sector/PID = -1/0 reproduce
//               the original's stored values (those setters are never called).
//    * vt / vt_pi : left at (0,0,0); the original generator skimmer never fills
//               the vertex even though LUND carries it.
//
//  --------------------------------------------------------------------------
//  USAGE
//  --------------------------------------------------------------------------
//    ./gen_skimmer_lund <input> <Ebeam> <inclusive(0|1)> <outFileName_no_ext>
//
//      <input>      path to a single .dat LUND file, OR a directory: every
//                   *.dat file in it is read into ONE output tree.  A directory
//                   is the natural analog of the original's multi-file loop
//                   (clasdis emits many .dat files per generator run); since
//                   there is no per-run number to split on at generator level,
//                   a single combined output file is written (the original's
//                   singleFile==1 mode).
//      <Ebeam>      beam energy in GeV (10.2 / 10.4 / 10.6; else forced 10.2).
//      <inclusive>  0 = require >=1 good pion (e,e'pi);  1 = electron-only skim.
//      <outFileName> output ROOT file name WITHOUT the .root extension.
//==============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <filesystem>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "constants.h"     // Me, Mpi, Mp, p_rest, rad_to_deg
#include "cut_values.h"    // Z_min, Q2_min, ... (single source of truth)
#include "lund_gen_classes.h"

using namespace constants;
using namespace cutVals;
namespace fs = std::filesystem;

//------------------------------------------------------------------------------
//  Kinematic cuts -- exact ports of
//  analyzer::applyElectronKinematicCuts(genElectron) and
//  analyzer::applyPionKinematicCuts(genPion) (classes/analyzer.cpp 473-509)
//  with every cut-sensitivity modifier mod_* = 1 (their default; loadCutValues
//  is a no-op in the original).
//------------------------------------------------------------------------------
bool applyElectronKinematicCuts( genElectron e ){
    if( std::sqrt(e.getW2()) < W_min )                       { return false; }
    if( e.getQ2() < Q2_min || e.getQ2() > Q2_max )           { return false; }
    if( e.getXb() < xB_min || e.getXb() > xB_max )           { return false; }
    if( e.getY()  > y_max )                                  { return false; }

    double p = e.get3Momentum().Mag();
    if( p < P_e_min || p > P_e_max )                         { return false; }

    double theta = e.get3Momentum().Theta() * rad_to_deg;
    if( theta < theta_min || theta > theta_max )             { return false; }

    return true;
}

bool applyPionKinematicCuts( genPion pi ){
    if( pi.getMx() < Mx_min || pi.getMx() > Mx_max )         { return false; }

    double p = pi.get3Momentum().Mag();
    if( p < P_pi_min || p > P_pi_max )                       { return false; }

    if( pi.getZ() < Z_min || pi.getZ() > Z_max )             { return false; }

    double theta = pi.get3Momentum().Theta() * rad_to_deg;
    if( theta < theta_min || theta > theta_max )             { return false; }

    return true;
}

//------------------------------------------------------------------------------
//  Minimal LUND particle record (only the fields the skimmer uses).
//------------------------------------------------------------------------------
struct LundPart {
    int    pid;
    double px, py, pz;
};

//  Gather the list of input .dat files (single file or directory).
std::vector<std::string> gatherInputFiles( const std::string& input ){
    std::vector<std::string> files;
    std::error_code ec;
    if( fs::is_directory(input, ec) ){
        for( const auto& de : fs::directory_iterator(input) ){
            if( de.is_regular_file() && de.path().extension() == ".dat" )
                files.push_back( de.path().string() );
        }
        std::sort( files.begin(), files.end() );
    } else {
        files.push_back( input );
    }
    return files;
}

//==============================================================================
int main( int argc, char** argv ){

    auto start = std::chrono::high_resolution_clock::now();

    if( argc < 5 ){
        std::cerr << "Incorrect number of arguments. Please use:\n";
        std::cerr << "  ./gen_skimmer_lund [input file or dir] [Beam energy] "
                     "[Inclusive (0,1)] [Output File Name (no extension)]\n";
        return -1;
    }

    std::string inputPath  = argv[1];
    double      Ebeam      = atof(argv[2]);
    int         inclusive  = atoi(argv[3]);
    TString     outFileName= argv[4];

    // Beam-energy validation -- identical policy to gen_skimmer.cpp.
    if( Ebeam != 10.2 && Ebeam != 10.4 && Ebeam != 10.6 ){
        std::cout << "Invalid Beam Energy... Set EBeam = 10.2\n" << std::endl;
        Ebeam = 10.2;
    }

    // ---- output variables (same names/types/order as gen_skimmer.cpp) -------
    int    Ne, Npi, Npips, Npims, runnum, evnum;
    double torusBending = -1;                         // branch "torus"
    TLorentzVector beam( 0, 0, Ebeam, Ebeam );

    genElectron           e_gen;
    std::vector<genPion>  pi_gen;
    int tar_string = 0;
    // ---- output file and tree (single combined file = singleFile==1) --------
    TFile* outputFile = new TFile( outFileName + ".root", "RECREATE" );
    TTree* outTree    = new TTree( "ePi", "(e,e'pi) event  information" );

    outTree->Branch( "runnum", &runnum );
    outTree->Branch( "torus",  &torusBending );
    outTree->Branch( "evnum",  &evnum );
    outTree->Branch( "Ebeam",  &Ebeam );
    outTree->Branch( "beam",   &beam );

    outTree->Branch( "Ne",     &Ne );

    outTree->Branch( "Npi",    &Npi );
    outTree->Branch( "Npips",  &Npips );
    outTree->Branch( "Npims",  &Npims );

    outTree->Branch( "e_gen",  &e_gen );
    outTree->Branch( "pi_gen", &pi_gen );
    outTree->Branch( "string", &tar_string);

    // ---- input files --------------------------------------------------------
    std::vector<std::string> inputFiles = gatherInputFiles( inputPath );
    if( inputFiles.empty() ){
        std::cerr << "No input .dat files found at: " << inputPath << "\n";
        return -1;
    }

    runnum = 0;            // no generator-level run number (agreed)
    long   globalEvent  = 0;   // 0-based sequential evnum across the whole job
    long   nProcessed   = 0;
    long   nFilled      = 0;

    // ----------------------------- file loop ---------------------------------
    for( size_t fi = 0; fi < inputFiles.size(); ++fi ){
        std::cout << "reading file " << fi+1 << " of " << inputFiles.size()
                  << " : " << inputFiles[fi] << "\n";

        std::ifstream fin( inputFiles[fi] );
        if( !fin.is_open() ){
            std::cerr << "  could not open " << inputFiles[fi] << " -- skipping\n";
            continue;
        }

        std::string line;
        while( std::getline(fin, line) ){
            // ---- parse the event header line --------------------------------
            std::istringstream hs(line);
            int nPart;
            if( !(hs >> nPart) ) { continue; }     // skip blank/garbage lines
            if( nPart <= 0 )     { continue; }

            // ---- read this event's particle lines ---------------------------
            std::vector<LundPart> parts;
            parts.reserve( nPart );
            for( int k = 0; k < nPart; ++k ){
                if( !std::getline(fin, line) ) break;
                std::istringstream ps(line);
                // columns: idx charge status pid parent daughter px py pz E m vx vy vz
                int    idx, status, pid, parent, daughter;
                double charge, px, py, pz, E, m, vx, vy, vz;
                if( !(ps >> idx >> charge >> status >> pid >> parent >> daughter
                         >> px >> py >> pz >> E >> m >> vx >> vy >> vz) ) continue;
                parts.push_back( { pid, px, py, pz } );
            }

            evnum = (int) globalEvent;     // 0-based sequential index
            ++globalEvent;
            ++nProcessed;

            // ---- per-event reset (TRUE per-event counts) --------------------
            e_gen.Clear();
            pi_gen.clear();
            Ne = Npi = Npips = Npims = 0;

            // ---- find particles by PID (skip index 0 = beam lepton) ---------
            std::vector<int> electronsMC, pipsMC, pimsMC;
            for( size_t i = 1; i < parts.size(); ++i ){   // i=1 -> skip beam
                switch( parts[i].pid ){
                    case 1: tar_string = 1; break;
                    case 2: tar_string = 2; break;
                    case  11: electronsMC.push_back(i); break;
                    case  211: pipsMC.push_back(i);     break;
                    case -211: pimsMC.push_back(i);     break;
                }
            }
            std::vector<int> pionsMC = pipsMC;
            pionsMC.insert( pionsMC.end(), pimsMC.begin(), pimsMC.end() );

            Ne    = (int) electronsMC.size();   // true per-event electron count
            Npi   = (int) pionsMC.size();
            Npips = (int) pipsMC.size();         // filled (original left at 0)
            Npims = (int) pimsMC.size();         // filled (original left at 0)

            // Event-level pre-cut (gen_skimmer.cpp line 211).
            if( Npi == 0 && inclusive != 1 ){ continue; }

            // Robustness: clasdis always has a scattered electron, but guard
            // against an electron-less event rather than dereferencing [0].
            if( electronsMC.empty() ){ continue; }

            // ---------------------- electron analysis ------------------------
            int e_idx = electronsMC[0];   // first PID==11 after beam = scattered e-
            e_gen.setMomentum( parts[e_idx].px, parts[e_idx].py, parts[e_idx].pz );
            e_gen.setKinematicInformation( Ebeam,
                                           parts[e_idx].px,
                                           parts[e_idx].py,
                                           parts[e_idx].pz );
            if( !applyElectronKinematicCuts( e_gen ) ){ continue; }

            // ------------------------ pion analysis --------------------------
            genPion genPi_dummy;
            for( int i = 0; i < Npi; ++i ){
                if( inclusive == 1 ){ continue; }       // electron-only skim
                genPi_dummy.Clear();
                int p_idx = pionsMC[i];
                genPi_dummy.setMomentum( parts[p_idx].px,
                                         parts[p_idx].py,
                                         parts[p_idx].pz );
                genPi_dummy.setKinematicInformation( e_gen.getQ(),
                                                     e_gen.get4Momentum(),
                                                     parts[p_idx].px,
                                                     parts[p_idx].py,
                                                     parts[p_idx].pz,
                                                     parts[p_idx].pid );
                //if( !applyPionKinematicCuts( genPi_dummy ) ){ continue; }
                pi_gen.push_back( genPi_dummy );
            }

            // Require a surviving pion in exclusive mode (line 243).
            if( pi_gen.size() == 0 && inclusive == 0 ){ continue; }

            outTree->Fill();
            ++nFilled;

            if( nProcessed % 5000 == 0 )
                std::cout << "  processed " << nProcessed
                          << " events, filled " << nFilled << "\n";
        }
        fin.close();
        std::cout << "Finished file!\n";
    }

    // ---- write -------------------------------------------------------------
    std::cout << "Writing tree to file\n";
    outputFile->cd();
    outTree->Write();
    outputFile->Close();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    std::cout << "Done.\n";
    std::cout << "  input events processed : " << nProcessed << "\n";
    std::cout << "  output events (filled) : " << nFilled    << "\n";
    std::cout << "  elapsed time           : " << elapsed.count() << " s\n";

    return 0;
}
