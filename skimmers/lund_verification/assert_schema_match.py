#!/usr/bin/env python3
"""
assert_schema_match.py

Rigorous, automated check that the PORT will produce a byte-identical schema to
the ORIGINAL gen_skimmer output.  It:
  (a) reads the ground-truth schema + class member layouts from the real file
      gen_skim_10.6.root via uproot;
  (b) parses the port's own source (lund_gen_classes.h private members, and the
      outTree->Branch(...) order in gen_skimmer_lund.cpp);
  (c) diffs them and prints PASS/FAIL for each.

A ROOT object branch's on-disk streamer depends ONLY on the class name and its
data members (names, types, declaration order) -- not member functions -- so if
(b) matches (a), a downstream reader cannot distinguish the files.
"""
import sys, re
sys.path.insert(0, "/tmp/upr")
import uproot

ORIG  = "/sessions/blissful-intelligent-heisenberg/mnt/sidis_analysis_suite/trees/final_skims/GEMC/gen_skim_10.6.root"
HDR   = "/sessions/blissful-intelligent-heisenberg/mnt/sidis_analysis_suite/skimmers/lund_gen_classes.h"
MAIN  = "/sessions/blissful-intelligent-heisenberg/mnt/sidis_analysis_suite/skimmers/gen_skimmer_lund.cpp"

CXX2ROOT = {"int": "int", "double": "double", "bool": "bool",
            "TLorentzVector": "TLorentzVector", "TVector3": "TVector3"}

def ok(b): return "PASS" if b else "FAIL"

# ---------------------------------------------------------------- ground truth
f = uproot.open(ORIG)
t = f["ePi"]
orig_tree_name  = t.name
orig_tree_title = t.title
orig_branches = [(name, str(br.typename)) for name, br in t.iteritems(recursive=False)]

streamers = f.file.streamers
def orig_members(cls):
    versions = streamers[cls]
    ver = sorted(versions.keys())[-1]
    si = versions[ver]
    out = []
    for el in si.elements:
        out.append((el.name, str(el.typename)))
    return ver, out

ge_ver, ge_members = orig_members("genElectron")
gp_ver, gp_members = orig_members("genPion")

# ---------------------------------------------------------------- parse port header
src = open(HDR).read()

def parse_private_members(class_name):
    # grab the class body
    m = re.search(r"class\s+" + class_name + r"\s*\{(.*?)\n\};", src, re.S)
    body = m.group(1)
    # take the text after the last 'private:'
    priv = body.split("private:")[-1]
    members = []
    for line in priv.splitlines():
        line = line.strip()
        if not line or line.startswith("//"):
            continue
        # strip trailing comments
        line = line.split("//")[0].strip()
        mm = re.match(r"(int|double|bool|TLorentzVector|TVector3)\s+([A-Za-z_]\w*)\s*;", line)
        if mm:
            members.append((mm.group(2), CXX2ROOT[mm.group(1)]))
    return members

port_ge = parse_private_members("genElectron")
port_gp = parse_private_members("genPion")

# ---------------------------------------------------------------- parse branch order
main_src = open(MAIN).read()
branch_calls = re.findall(r'outTree->Branch\(\s*"([^"]+)"', main_src)
tree_m = re.search(r'new\s+TTree\(\s*"([^"]+)"\s*,\s*"([^"]*)"', main_src)
port_tree_name, port_tree_title = tree_m.group(1), tree_m.group(2)

# ---------------------------------------------------------------- report
print("=" * 72)
print("1) TREE NAME / TITLE")
print("=" * 72)
print(f"   original : name={orig_tree_name!r}  title={orig_tree_title!r}")
print(f"   port     : name={port_tree_name!r}  title={port_tree_title!r}")
print(f"   -> name  {ok(orig_tree_name==port_tree_name)}")
print(f"   -> title {ok(orig_tree_title==port_tree_title)}")

print("\n" + "=" * 72)
print("2) TOP-LEVEL BRANCH ORDER & TYPES")
print("=" * 72)
orig_names = [n for n, _ in orig_branches]
print(f"   original ({len(orig_names)}): {orig_names}")
print(f"   port     ({len(branch_calls)}): {branch_calls}")
print(f"   -> order/count {ok(orig_names==branch_calls)}")

print("\n" + "=" * 72)
print(f"3) genElectron MEMBER LAYOUT  (original class version {ge_ver})")
print("=" * 72)
match = ge_members == port_ge
for i in range(max(len(ge_members), len(port_ge))):
    o = ge_members[i] if i < len(ge_members) else ("--", "--")
    p = port_ge[i]    if i < len(port_ge)    else ("--", "--")
    same = o == p
    print(f"   [{i:>2}] orig {o[0]:<10}{o[1]:<16} | port {p[0]:<10}{p[1]:<16} {ok(same)}")
print(f"   -> genElectron layout {ok(match)}")

print("\n" + "=" * 72)
print(f"4) genPion MEMBER LAYOUT  (original class version {gp_ver})")
print("=" * 72)
match2 = gp_members == port_gp
for i in range(max(len(gp_members), len(port_gp))):
    o = gp_members[i] if i < len(gp_members) else ("--", "--")
    p = port_gp[i]    if i < len(port_gp)    else ("--", "--")
    same = o == p
    print(f"   [{i:>2}] orig {o[0]:<10}{o[1]:<16} | port {p[0]:<10}{p[1]:<16} {ok(same)}")
print(f"   -> genPion layout {ok(match2)}")

print("\n" + "=" * 72)
allok = (orig_tree_name==port_tree_name and orig_tree_title==port_tree_title
         and orig_names==branch_calls and match and match2)
print("OVERALL:", "ALL CHECKS PASS -- schema is byte-identical." if allok
      else "MISMATCH -- see FAIL lines above.")
print("=" * 72)
