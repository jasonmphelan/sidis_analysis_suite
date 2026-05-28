#!/bin/bash
# =============================================================================
# Theta_pi acceptance matching cut scan
# Scans matchCut2D_XX.root files (XX = 5, 10, 15, 20)
# For each: recompile, run corrections, run analysis, extract R_d
# =============================================================================

set -e  # exit on error

BASEDIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILDDIR="${BASEDIR}/build"
SCANDIR="${BASEDIR}/theta_scan"
CORRDIR="${BASEDIR}/data/correctionFiles"
HISTDIR="${SCANDIR}"

LOGFILE="${SCANDIR}/run_log.txt"
SUMMARY="${SCANDIR}/summary.txt"

# Input files for corrections (10.2 GeV combined)
REC_FILE="${BASEDIR}/trees/final_skims/GEMC/final_skim_10.2_rga.root"
GEN_FILE="${BASEDIR}/trees/final_skims/GEMC/gen_skim_10.2_rga.root"

# Correction parameters
MATCH_TYPE=1    # 2D acceptance matching
TARGET=0        # RGB/deuterium
THETA_CUT=0     # flat theta cut (inactive)

# Analysis parameters
BEAM=10.2
APPLY_CORR=3    # MC correction (acc + bin migration)
MAP=1
ACC_FILE="acceptanceMap_allE_final.root"
VAR_NAME="Mx"
BINS_VAR=1
VAR_MIN=1.7
VAR_MAX=3.7
BINNED_OUT=1
MX2_CUT=1.25
RHO_MODE=1

# Thresholds to scan
THRESHOLDS=(1 5 10 15 20)

# Initialize log and summary
echo "Theta_pi matching cut scan - $(date)" > "$LOGFILE"
echo "============================================" >> "$LOGFILE"
echo "# threshold  mean_Rd(z=0.3-0.7)  correction_file  ratio_file" > "$SUMMARY"

# -----------------------------------------------------------------------------
# Step 0: Recompile calcCorrections (only needed once after source edit)
# -----------------------------------------------------------------------------
echo "=== Recompiling ===" | tee -a "$LOGFILE"
cd "$BUILDDIR"
cmake .. >> "$LOGFILE" 2>&1
make calcCorrections 2>&1 | tee -a "$LOGFILE"
make makeRatio3D 2>&1 | tee -a "$LOGFILE"
echo "Compilation complete." | tee -a "$LOGFILE"
echo "" >> "$LOGFILE"

# -----------------------------------------------------------------------------
# Step 1: Run scan
# -----------------------------------------------------------------------------
for THR in "${THRESHOLDS[@]}"; do
    MATCH_FILE="matchCut2D_${THR}.root"
    CORR_OUT="${CORRDIR}/corrections_10.2_theta_scan_${THR}.root"
    RATIO_OUT="${HISTDIR}/ratio_Rd_theta_${THR}.root"
    RD_TXT="${HISTDIR}/Rd_theta_${THR}.txt"

    echo "" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"
    echo "=== Threshold: ${THR}% (${MATCH_FILE}) ===" | tee -a "$LOGFILE"
    echo "============================================" | tee -a "$LOGFILE"

    # --- Run corrections ---
    echo "--- Running calcCorrections ---" | tee -a "$LOGFILE"
    CORR_CMD="${BUILDDIR}/corrections/calcCorrections ${REC_FILE} ${GEN_FILE} ${CORR_OUT} ${MATCH_TYPE} ${TARGET} ${THETA_CUT} ${MATCH_FILE}"
    echo "CMD: ${CORR_CMD}" >> "$LOGFILE"
    ${CORR_CMD} 2>&1 | tee -a "$LOGFILE"
    CORR_EXIT=$?
    if [ $CORR_EXIT -ne 0 ] && [ $CORR_EXIT -ne 1 ]; then
        echo "ERROR: calcCorrections failed with exit code ${CORR_EXIT}" | tee -a "$LOGFILE"
        echo "Stopping scan." | tee -a "$LOGFILE"
        exit 1
    fi
    echo "Corrections written to: ${CORR_OUT}" | tee -a "$LOGFILE"

    # --- Run analysis ---
    echo "--- Running makeRatio3D ---" | tee -a "$LOGFILE"
    CORR_NAME="corrections_10.2_theta_scan_${THR}.root"
    RATIO_CMD="${BUILDDIR}/analysis/makeRatio3D ${BEAM} ${RATIO_OUT} ${MATCH_TYPE} ${APPLY_CORR} ${MAP} ${ACC_FILE} ${VAR_NAME} ${BINS_VAR} ${VAR_MIN} ${VAR_MAX} ${BINNED_OUT} ${CORR_NAME} ${MX2_CUT} ${RHO_MODE} ${TARGET} ${MATCH_FILE}"
    echo "CMD: ${RATIO_CMD}" >> "$LOGFILE"
    ${RATIO_CMD} 2>&1 | tee -a "$LOGFILE"
    RATIO_EXIT=$?
    if [ $RATIO_EXIT -ne 0 ] && [ $RATIO_EXIT -ne 1 ]; then
        echo "ERROR: makeRatio3D failed with exit code ${RATIO_EXIT}" | tee -a "$LOGFILE"
        echo "Stopping scan." | tee -a "$LOGFILE"
        exit 1
    fi
    echo "Ratios written to: ${RATIO_OUT}" | tee -a "$LOGFILE"

    # --- Extract R_d values ---
    echo "--- Extracting R_d ---" | tee -a "$LOGFILE"
    root -l -b -q "${SCANDIR}/extract_Rd.C(\"${RATIO_OUT}\", \"${RD_TXT}\")" 2>&1 | tee -a "$LOGFILE"

    # --- Compute mean R_d for z = 0.3-0.7 and append to summary ---
    if [ -f "$RD_TXT" ]; then
        # Average R_d over z bins with center <= 0.7 (first 8 bins: 0.325 to 0.675)
        MEAN_RD=$(awk 'NR>1 && $1 <= 0.7 { sum+=$2; n++ } END { if(n>0) printf "%.4f", sum/n; else print "NaN" }' "$RD_TXT")
        echo "${THR}  ${MEAN_RD}  ${CORR_OUT}  ${RATIO_OUT}" >> "$SUMMARY"
        echo "" | tee -a "$LOGFILE"
        echo "R_d values for threshold ${THR}%:" | tee -a "$LOGFILE"
        cat "$RD_TXT" | tee -a "$LOGFILE"
        echo "Mean R_d (z=0.3-0.7): ${MEAN_RD}" | tee -a "$LOGFILE"
    else
        echo "WARNING: ${RD_TXT} not found" | tee -a "$LOGFILE"
        echo "${THR}  NaN  ${CORR_OUT}  ${RATIO_OUT}" >> "$SUMMARY"
    fi

    echo "" | tee -a "$LOGFILE"
    echo "=== Threshold ${THR}% complete ===" | tee -a "$LOGFILE"
done

echo "" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "=== Scan complete ===" | tee -a "$LOGFILE"
echo "============================================" | tee -a "$LOGFILE"
echo "" | tee -a "$LOGFILE"
echo "Summary:" | tee -a "$LOGFILE"
cat "$SUMMARY" | tee -a "$LOGFILE"

# --- Print comparison table ---
echo "" | tee -a "$LOGFILE"
echo "=== Comparison Table ===" | tee -a "$LOGFILE"
root -l -b -q "${SCANDIR}/compare_Rd.C" 2>&1 | tee -a "$LOGFILE"
