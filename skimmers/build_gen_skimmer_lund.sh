#!/usr/bin/env bash
#==============================================================================
#  build_gen_skimmer_lund.sh
#
#  Standalone build for the LUND (clasdis) skimmer port.
#
#  This build is intentionally SELF-CONTAINED: it needs only ROOT (root-config
#  on your PATH) and a C++17 compiler.  It does NOT need clas12root, hipo, or
#  the project CMake build, because gen_skimmer_lund.cpp reads plain-text LUND
#  files and uses the standalone genElectron/genPion classes in
#  lund_gen_classes.h (no CLAS12 dependencies).
#
#  It does NOT touch the existing CMakeLists.txt or any pre-existing source.
#
#  Usage:
#     cd skimmers
#     ./build_gen_skimmer_lund.sh
#  Produces:
#     ./gen_skimmer_lund          (executable, in the skimmers directory)
#
#  Run example:
#     ./gen_skimmer_lund ../gemc/eventfiles/clasdisde.00.e10.200.emn0.75tmn.09.xs0000.dat \
#                        10.2 0 my_lund_skim
#==============================================================================
set -euo pipefail

# --- locate things relative to this script -----------------------------------
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
INCLUDE_DIR="$( cd "${HERE}/../include" && pwd )"

# --- sanity: ROOT must be available ------------------------------------------
if ! command -v root-config >/dev/null 2>&1 ; then
    echo "ERROR: root-config not found on PATH. Set up ROOT first (e.g. 'source <root>/bin/thisroot.sh')." >&2
    exit 1
fi

CXX="$(root-config --cxx)"
ROOT_CFLAGS="$(root-config --cflags)"
ROOT_LIBS="$(root-config --libs)"

echo "Using compiler : ${CXX}"
echo "ROOT cflags    : ${ROOT_CFLAGS}"
echo "include dir    : ${INCLUDE_DIR}"

cd "${HERE}"

# --- 1. generate the ROOT dictionary for the standalone classes --------------
#     -I flags must let rootcling find lund_gen_classes.h AND constants.h.
echo ">>> generating dictionary (rootcling)"
rootcling -f lund_gen_dict.cxx \
    -I"${HERE}" -I"${INCLUDE_DIR}" \
    lund_gen_classes.h lund_gen_classes_LinkDef.h

# --- 2. compile + link the executable ----------------------------------------
echo ">>> compiling gen_skimmer_lund"
${CXX} -std=c++17 -O2 ${ROOT_CFLAGS} \
    -I"${HERE}" -I"${INCLUDE_DIR}" \
    gen_skimmer_lund.cpp lund_gen_dict.cxx \
    ${ROOT_LIBS} -lEG \
    -o gen_skimmer_lund

echo ">>> done: ${HERE}/gen_skimmer_lund"
