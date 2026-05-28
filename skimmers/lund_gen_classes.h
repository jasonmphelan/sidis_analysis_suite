#ifndef LUND_GEN_CLASSES_H
#define LUND_GEN_CLASSES_H
//==============================================================================
//  lund_gen_classes.h
//
//  Standalone re-implementations of genElectron and genPion for the clasdis
//  LUND skimmer (gen_skimmer_lund.cpp).
//
//  WHY THIS FILE EXISTS
//  --------------------
//  The output tree of gen_skimmer.cpp stores two object branches,
//      e_gen  : genElectron
//      pi_gen : vector<genPion>
//  whose on-disk ROOT StreamerInfo must be IDENTICAL to the original so that a
//  downstream analysis program cannot tell which skimmer wrote the file.  A
//  ROOT object branch's streamer is determined ONLY by the class name plus its
//  data members (names, types, declaration order) -- NOT by its member
//  functions.  The data members below are therefore copied verbatim, in the
//  same order, from include/genElectron.h and include/genPion.h.  Verified
//  against real gen_skimmer output (trees/final_skims/GEMC/gen_skim_10.6.root):
//  both classes are stored as class version 1 with exactly these members.
//
//  The ONLY behavioural change versus the originals is the momentum setter:
//  the original setMomentum(clas12::mcpar_ptr) pulls in the entire CLAS12
//  framework, which is meaningless for a plain-text LUND reader.  It is
//  replaced here by a direct (px,py,pz) setter.  All kinematic formulas are
//  copied verbatim from classes/genElectron.cpp and classes/genPion.cpp.
//
//  Member functions differ from the originals (extra/changed setters) but that
//  has NO effect on the ROOT streamer, hence none on tree compatibility.
//==============================================================================

#include "TLorentzVector.h"
#include "TVector3.h"
#include "constants.h"   // constants:: Me, Mpi, Mp, p_rest, rad_to_deg (project header, no CLAS12 deps)
#include <cmath>

//------------------------------------------------------------------------------
//  genElectron  (data-member-compatible with include/genElectron.h)
//------------------------------------------------------------------------------
class genElectron {
    public:
        genElectron(){ Clear(); }
        ~genElectron(){}

        void Clear(){
            Sector    = -1;
            PID       = 0;
            Charge    = 0;

            Q2        = 0;
            Omega     = 0;
            Xb        = 0;
            W2        = 0;
            y         = 0;

            Selection = false;

            e4 = TLorentzVector();
            e3 = TVector3();
            q  = TLorentzVector();
            vt = TVector3();
        }

        // ---- particle ID accessors (same API as original) -------------------
        void setSector (int v){ Sector = v; }
        void setPID    (int v){ PID    = v; }
        void setCharge (int v){ Charge = v; }
        int  getSector (void){ return Sector; }
        int  getPID    (void){ return PID;    }
        int  getCharge (void){ return Charge; }

        // ---- kinematics accessors -------------------------------------------
        void setQ  (TLorentzVector iQ){ q = iQ; }
        TLorentzVector get4Momentum (void){ return e4; }
        TVector3       get3Momentum (void){ return e3; }
        TLorentzVector getQ         (void){ return q;  }

        void setQ2    (double v){ Q2    = v; }
        void setOmega (double v){ Omega = v; }
        void setXb    (double v){ Xb    = v; }
        void setW2    (double v){ W2    = v; }
        void setY     (double v){ y     = v; }
        double getQ2    (void){ return Q2;    }
        double getOmega (void){ return Omega; }
        double getXb    (void){ return Xb;    }
        double getW2    (void){ return W2;    }
        double getY     (void){ return y;     }

        void setSelection (bool v){ Selection = v; }
        bool getSelection (void){ return Selection; }

        // ---- LUND-port momentum setter (replaces mcpar_ptr version) ----------
        //  Mirrors genElectron::setMomentum: build e4 from px,py,pz with the
        //  fixed electron mass and take e3 as its 3-vector.
        void setMomentum (double px, double py, double pz){
            e4.SetXYZM( px, py, pz, constants::Me );
            e3 = e4.Vect();
        }

        // ---- verbatim port of genElectron::setKinematicInformation -----------
        void setKinematicInformation (double Ebeam, double px, double py, double pz){
            TLorentzVector beam( 0, 0, Ebeam, Ebeam );
            setMomentum( px, py, pz );

            setQ   ( beam - e4 );
            setQ2  ( (double) -q.Mag2() );
            setOmega( q.E() );
            setXb  ( Q2 / (2. * constants::Mp * q.E()) );
            setY   ( Omega / Ebeam );
            setW2  ( (double)((constants::p_rest + q).Mag2()) );
        }

    private:
        int Sector;
        int PID;
        int Charge;

        TLorentzVector e4;
        TVector3       e3;
        TLorentzVector q;
        TVector3       vt;

        double Q2;
        double Omega;
        double Xb;
        double W2;
        double y;

        bool Selection;
};

//------------------------------------------------------------------------------
//  genPion  (data-member-compatible with include/genPion.h)
//------------------------------------------------------------------------------
class genPion {
    public:
        genPion(){ Clear(); }
        ~genPion(){}

        //  NOTE: the original genPion::Clear() does NOT reset Sector or Charge,
        //  leaving Sector indeterminate in the original output.  Per the agreed
        //  design we pin Sector to a deterministic value (-1, matching
        //  genElectron's "no sector" convention) instead of emitting
        //  nondeterministic stack garbage.  PID is left at 0 exactly as the
        //  original (setPID is never called for pions).
        void Clear(){
            Sector    = -1;   // deterministic (original: uninitialised)
            PID       = 0;
            Charge    = 0;

            Z         = 0;
            Mx        = 0;
            xF        = 0;
            eta       = 0;

            Selection = false;

            pi4   = TLorentzVector();
            pi3   = TVector3();
            pi_q  = TLorentzVector();
            vt_pi = TVector3();
        }

        // ---- particle ID accessors ------------------------------------------
        void setSector (int v){ Sector = v; }
        void setPID    (int v){ PID    = v; }
        void setCharge (int v){ Charge = v; }
        int  getSector (void){ return Sector; }
        int  getPID    (void){ return PID;    }
        int  getCharge (void){ return Charge; }

        // ---- kinematics accessors -------------------------------------------
        TLorentzVector get4Momentum (void){ return pi4;  }
        TVector3       get3Momentum (void){ return pi3;  }
        TLorentzVector getPi_q      (void){ return pi_q; }
        TVector3       getVt        (void){ return vt_pi;}

        void setZ   (double v){ Z   = v; }
        void setMx  (double v){ Mx  = v; }
        void setXf  (double v){ xF  = v; }
        void setEta (double v){ eta = v; }
        double getZ   (void){ return Z;   }
        double getMx  (void){ return Mx;  }
        double getXf  (void){ return xF;  }
        double getEta (void){ return eta; }

        void setSelection (bool v){ Selection = v; }
        bool getSelection (void){ return Selection; }

        // ---- LUND-port momentum setter (replaces mcpar_ptr version) ----------
        void setMomentum (double px, double py, double pz){
            pi4.SetXYZM( px, py, pz, constants::Mpi );
            pi3 = pi4.Vect();
        }

        // ---- verbatim port of genPion::setPi_q -------------------------------
        //  Boost/rotate into the q-Pe system: q is the z axis, Pe in x-z plane.
        void setPi_q (TLorentzVector q, TLorentzVector pe){
            TVector3 pi3_temp = pi3;
            pi3_temp.RotateZ( -q.Phi()   );
            pi3_temp.RotateY( -q.Theta() );
            pi3_temp.RotateZ( -pe.Phi()  );
            pi_q.SetVectM( pi3_temp, constants::Mpi );
        }

        // ---- verbatim port of genPion::setKinematicInformation ---------------
        //  pid is needed for Charge = (int)(pid/211); it came from rp->getPid()
        //  in the original.  xF is intentionally left at 0 (the original line
        //  that would set it is commented out).
        void setKinematicInformation (TLorentzVector q, TLorentzVector pe,
                                      double px, double py, double pz, int pid){
            if( pi4.P() == 0 ){ setMomentum( px, py, pz ); }

            setPi_q( q, pe );

            setCharge( (int) (pid / 211) );

            setZ ( pi4.E() / q.E() );

            setMx( ( q + constants::p_rest - pi4 ).Mag() );

            setEta( 0.5 * log( (pi_q.E() + pi_q.Pz()) /
                               (pi_q.E() - pi_q.Pz()) ) );
        }

    private:
        int Sector;
        int PID;
        int Charge;

        TLorentzVector pi4;
        TVector3       pi3;
        TLorentzVector pi_q;
        TVector3       vt_pi;

        double Z;
        double Mx;
        double eta;
        double xF;

        bool Selection;
};

#endif // LUND_GEN_CLASSES_H
