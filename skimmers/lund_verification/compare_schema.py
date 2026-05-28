#!/usr/bin/env python3
"""Extract the exact on-disk schema + streamer layout of the ORIGINAL
gen_skimmer output (gen_skim_10.6.root) and print it as a TTree::Print-style
report, so it can be compared branch-by-branch with what gen_skimmer_lund.cpp
will write."""
import sys
sys.path.insert(0, "/tmp/upr")
import uproot

ORIG = "/sessions/blissful-intelligent-heisenberg/mnt/sidis_analysis_suite/trees/final_skims/GEMC/gen_skim_10.6.root"

f = uproot.open(ORIG)
print("FILE:", ORIG)
print("keys:", f.keys(), "\n")

t = f["ePi"]
print("=" * 72)
print(f"TTree '{t.name}'   title = {t.title!r}")
print(f"  entries = {t.num_entries}")
print("=" * 72)
print(f"{'#':>2}  {'branch':<10} {'typename':<28} {'interpretation'}")
print("-" * 72)
for i, (name, br) in enumerate(t.iteritems()):
    print(f"{i:>2}  {name:<10} {str(br.typename):<28} {br.interpretation}")

# ---- streamer info for the two object classes ----
print("\n" + "=" * 72)
print("STREAMERINFO for genElectron / genPion (class version + members)")
print("=" * 72)
streamers = f.file.streamers
for cls in ("genElectron", "genPion"):
    if cls not in streamers:
        print(f"  {cls}: NOT FOUND in streamers")
        continue
    versions = streamers[cls]
    for ver, si in versions.items():
        print(f"\n  class {cls}  (version {ver}):")
        for el in si.member_names:
            print(f"      {el}")
        # full element typing
        try:
            for el in si.elements:
                print(f"        - {el.name:<10} type={el.typename}")
        except Exception:
            pass
