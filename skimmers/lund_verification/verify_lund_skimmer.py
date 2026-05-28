#!/usr/bin/env python3
"""
verify_lund_skimmer.py  --  Step-4 verification harness for gen_skimmer_lund.cpp

The sandbox has no ROOT / clas12root, so we cannot build the C++ skimmer here.
Instead this script re-implements, in pure numpy, the EXACT same parse +
kinematics + cuts that gen_skimmer_lund.cpp performs, and runs it on a real
clasdis LUND file.  It reports:

  * input event count and surviving (filled) event count
  * a particle-record audit proving the "skip index-0, select by PID 11/211/-211"
    rule is equivalent to selecting only final-state e-/pi+/pi- (i.e. no
    intermediate PYTHIA particle ever carries PID 11/211/-211)
  * a full spot-check of one surviving event: every cut variable + threshold

All TLorentzVector / TVector3 operations are replicated to match ROOT exactly.
"""

import sys, math

# ---------------------------------------------------------------- constants
Me  = 0.00051099895
Mpi = 0.139570
Mp  = 0.938272
RAD2DEG = 180.0 / math.pi

# ---------------------------------------------------------------- cut values
Z_min, Z_max   = 0.3, 1.0
Q2_min, Q2_max = 2.0, 8.0
xB_min, xB_max = 0.1, 0.66
y_max          = 0.75
W_min          = 2.5
theta_min, theta_max = 5.0, 35.0
P_pi_min, P_pi_max   = 1.25, 5.0
P_e_min,  P_e_max    = 3.0, 10.6
Mx_min, Mx_max       = 1.7, 5.0

# ---------------------------------------------------------------- 4-vector helpers
def lv_from_pxpypzm(px, py, pz, m):
    """ROOT TLorentzVector::SetXYZM -> (px,py,pz, E=sqrt(p^2+m^2))."""
    E = math.sqrt(px*px + py*py + pz*pz + m*m)
    return [px, py, pz, E]

def lv_sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]]

def lv_add(a, b):
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]]

def lv_mag2(v):
    """ROOT Mag2 = E^2 - px^2 - py^2 - pz^2."""
    return v[3]*v[3] - v[0]*v[0] - v[1]*v[1] - v[2]*v[2]

def lv_mag(v):
    m2 = lv_mag2(v)
    return math.sqrt(m2) if m2 >= 0 else -math.sqrt(-m2)

def v3_mag(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def v3_theta(v):
    return math.atan2(math.sqrt(v[0]*v[0] + v[1]*v[1]), v[2])

def v3_phi(v):
    return math.atan2(v[1], v[0])

def lv_phi(v):
    return math.atan2(v[1], v[0])

def lv_theta(v):
    return math.atan2(math.sqrt(v[0]*v[0] + v[1]*v[1]), v[2])

def rotate_z(vec, a):
    s, c = math.sin(a), math.cos(a)
    x, y, z = vec
    return [c*x - s*y, s*x + c*y, z]

def rotate_y(vec, a):
    s, c = math.sin(a), math.cos(a)
    x, y, z = vec
    return [c*x + s*z, y, -s*x + c*z]

# ---------------------------------------------------------------- kinematics
def electron_kin(Ebeam, px, py, pz):
    e4   = lv_from_pxpypzm(px, py, pz, Me)
    beam = [0.0, 0.0, Ebeam, Ebeam]
    q    = lv_sub(beam, e4)
    Q2   = -lv_mag2(q)
    nu   = q[3]
    Xb   = Q2 / (2.0 * Mp * nu)
    y    = nu / Ebeam
    p_rest = [0.0, 0.0, 0.0, Mp]
    W2   = lv_mag2(lv_add(p_rest, q))
    e3   = [px, py, pz]
    return dict(e4=e4, e3=e3, q=q, Q2=Q2, nu=nu, Xb=Xb, y=y, W2=W2)

def electron_pass(k):
    if math.sqrt(k["W2"]) < W_min:                     return False
    if k["Q2"] < Q2_min or k["Q2"] > Q2_max:           return False
    if k["Xb"] < xB_min or k["Xb"] > xB_max:           return False
    if k["y"]  > y_max:                                return False
    p = v3_mag(k["e3"])
    if p < P_e_min or p > P_e_max:                     return False
    th = v3_theta(k["e3"]) * RAD2DEG
    if th < theta_min or th > theta_max:               return False
    return True

def pion_kin(q, pe, px, py, pz, pid):
    pi4  = lv_from_pxpypzm(px, py, pz, Mpi)
    pi3  = [px, py, pz]
    # setPi_q: rotate into q-frame
    t = list(pi3)
    t = rotate_z(t, -lv_phi(q))
    t = rotate_y(t, -lv_theta(q))
    t = rotate_z(t, -lv_phi(pe))
    Epq = math.sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2] + Mpi*Mpi)  # SetVectM E
    pi_q = [t[0], t[1], t[2], Epq]
    charge = int(pid / 211)
    Z  = pi4[3] / q[3]
    p_rest = [0.0, 0.0, 0.0, Mp]
    Mx = lv_mag(lv_sub(lv_add(q, p_rest), pi4))
    eta = 0.5 * math.log((pi_q[3] + pi_q[2]) / (pi_q[3] - pi_q[2]))
    return dict(pi4=pi4, pi3=pi3, pi_q=pi_q, charge=charge, Z=Z, Mx=Mx, eta=eta)

def pion_pass(k):
    if k["Mx"] < Mx_min or k["Mx"] > Mx_max:           return False
    p = v3_mag(k["pi3"])
    if p < P_pi_min or p > P_pi_max:                   return False
    if k["Z"] < Z_min or k["Z"] > Z_max:               return False
    th = v3_theta(k["pi3"]) * RAD2DEG
    if th < theta_min or th > theta_max:               return False
    return True

# ---------------------------------------------------------------- LUND parsing
def parse_events(path):
    """Yield (parts, raw_parts) per event. parts = list of (pid,px,py,pz);
    raw_parts keeps (idx,status,pid) for the audit."""
    with open(path) as f:
        it = iter(f)
        for line in it:
            toks = line.split()
            if not toks:
                continue
            try:
                npart = int(toks[0])
            except ValueError:
                continue
            if npart <= 0:
                continue
            parts, raw = [], []
            for _ in range(npart):
                try:
                    pl = next(it)
                except StopIteration:
                    break
                pt = pl.split()
                if len(pt) < 14:
                    continue
                idx    = int(pt[0])
                status = int(pt[2])
                pid    = int(pt[3])
                px, py, pz = float(pt[6]), float(pt[7]), float(pt[8])
                parts.append((pid, px, py, pz))
                raw.append((idx, status, pid))
            yield parts, raw

# ---------------------------------------------------------------- main
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "/sessions/blissful-intelligent-heisenberg/mnt/sidis_analysis_suite/gemc/eventfiles/clasdisde.00.e10.200.emn0.75tmn.09.xs0000.dat"
    Ebeam = float(sys.argv[2]) if len(sys.argv) > 2 else 10.2
    inclusive = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    print(f"Input file : {path}")
    print(f"Ebeam={Ebeam}  inclusive={inclusive}\n")

    n_proc = 0
    n_filled = 0

    # audit: which (status) values appear for selected PIDs among index>=1?
    sel_status = {11: {}, 211: {}, -211: {}}
    # also audit: are there PID 11/211/-211 at index 0 other than beam? (sanity)
    idx0_pids = {}

    spot = None  # first surviving event details

    for parts, raw in parse_events(path):
        n_proc += 1

        # --- audit particle records ---
        if raw:
            p0 = raw[0][2]
            idx0_pids[p0] = idx0_pids.get(p0, 0) + 1
        for j in range(1, len(parts)):
            pid = parts[j][0]
            if pid in sel_status:
                st = raw[j][1]
                sel_status[pid][st] = sel_status[pid].get(st, 0) + 1

        # --- find particles by PID, skip index 0 (beam) ---
        electrons, pips, pims = [], [], []
        for j in range(1, len(parts)):
            pid = parts[j][0]
            if   pid ==  11: electrons.append(j)
            elif pid ==  211: pips.append(j)
            elif pid == -211: pims.append(j)
        pions = pips + pims

        Ne, Npi, Npips, Npims = len(electrons), len(pions), len(pips), len(pims)

        if Npi == 0 and inclusive != 1:
            continue
        if not electrons:
            continue

        # --- electron ---
        e_idx = electrons[0]
        _, epx, epy, epz = parts[e_idx]
        ek = electron_kin(Ebeam, epx, epy, epz)
        if not electron_pass(ek):
            continue

        # --- pions ---
        good_pions = []
        for pidx in pions:
            if inclusive == 1:
                continue
            ppid, ppx, ppy, ppz = parts[pidx]
            pk = pion_kin(ek["q"], ek["e4"], ppx, ppy, ppz, ppid)
            if not pion_pass(pk):
                continue
            good_pions.append((pidx, ppid, pk))

        if len(good_pions) == 0 and inclusive == 0:
            continue

        n_filled += 1
        if spot is None:
            spot = (n_proc - 1, Ne, Npi, Npips, Npims, e_idx, ek, good_pions, parts)

    # ---------------------------------------------------------------- report
    print("=" * 70)
    print("PARTICLE-RECORD AUDIT")
    print("=" * 70)
    print("PID at index 0 (the skipped 'beam' slot), counts over all events:")
    for pid, c in sorted(idx0_pids.items()):
        print(f"    pid {pid:>6} : {c}")
    print("\nFor selected PIDs among index>=1, status -> count:")
    print("  (clasdis/PYTHIA status: 1 = final state, 11/21 = doc/intermediate)")
    for pid in (11, 211, -211):
        sts = sel_status[pid]
        pretty = ", ".join(f"status {s}: {n}" for s, n in sorted(sts.items()))
        print(f"    pid {pid:>5} : {pretty if pretty else '(none)'}")
    nonfinal = any(s != 1 for pid in (11, 211, -211) for s in sel_status[pid])
    print("\n  => " + ("WARNING: some selected particles are NOT status 1!"
                        if nonfinal else
                        "All PID 11/211/-211 found at index>=1 are status 1 "
                        "(final state). Skip-index-0 + PID select == final-state select. OK."))

    print()
    print("=" * 70)
    print("EVENT COUNTS")
    print("=" * 70)
    print(f"  input events processed : {n_proc}")
    print(f"  output events (filled) : {n_filled}")
    if n_proc:
        print(f"  acceptance             : {100.0*n_filled/n_proc:.3f} %")

    # ---------------------------------------------------------------- spot check
    if spot is None:
        print("\nNo surviving events to spot-check.")
        return
    (evnum, Ne, Npi, Npips, Npims, e_idx, ek, good_pions, parts) = spot
    print()
    print("=" * 70)
    print(f"SPOT-CHECK of first surviving event (evnum={evnum})")
    print("=" * 70)
    print(f"  Ne={Ne}  Npi={Npi}  Npips={Npips}  Npims={Npims}")
    print(f"  scattered e- is parts[{e_idx}] (0-based, after skipping beam)")
    print("\n  ELECTRON kinematics & cuts:")
    W = math.sqrt(ek["W2"])
    p_e = v3_mag(ek["e3"])
    th_e = v3_theta(ek["e3"]) * RAD2DEG
    def chk(name, val, lo, hi):
        ok = (lo is None or val >= lo) and (hi is None or val <= hi)
        rng = f"[{lo}, {hi}]" if lo is not None and hi is not None else \
              (f">= {lo}" if hi is None else f"<= {hi}")
        print(f"     {name:<10}= {val:>10.4f}   need {rng:<16} {'PASS' if ok else 'FAIL'}")
    chk("W",   W,        W_min, None)
    chk("Q2",  ek["Q2"], Q2_min, Q2_max)
    chk("xB",  ek["Xb"], xB_min, xB_max)
    chk("y",   ek["y"],  None, y_max)
    chk("|p_e|", p_e,    P_e_min, P_e_max)
    chk("theta_e", th_e, theta_min, theta_max)

    print("\n  PION(S) kinematics & cuts:")
    for (pidx, ppid, pk) in good_pions:
        p_pi = v3_mag(pk["pi3"])
        th_pi = v3_theta(pk["pi3"]) * RAD2DEG
        print(f"   - parts[{pidx}] pid={ppid} charge={pk['charge']}:")
        chk("Mx",      pk["Mx"], Mx_min, Mx_max)
        chk("|p_pi|",  p_pi,     P_pi_min, P_pi_max)
        chk("Z",       pk["Z"],  Z_min, Z_max)
        chk("theta_pi",th_pi,    theta_min, theta_max)
        print(f"     eta       = {pk['eta']:>10.4f}   (stored, no cut)")

if __name__ == "__main__":
    main()
