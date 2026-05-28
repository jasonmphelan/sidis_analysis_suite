#!/bin/bash
# merge_gemc_data.sh
# Merges a GEMC simulation skim with a data run skim using hadd.
#
# Usage:
#   ./merge_gemc_data.sh [output_file] [gemc_file] [data_file]
#
# Example:
#   ./merge_gemc_data.sh merged.root \
#       ../trees/final_skims/GEMC/final_skim_10.2.root \
#       ../trees/final_skims/10.2/final_skim_10637.root

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [output_file] [gemc_file] [data_file]"
    exit 1
fi

OUT=$1
GEMC=$2
DATA=$3

if [ ! -f "$GEMC" ]; then
    echo "ERROR: GEMC file not found: $GEMC"
    exit 1
fi
if [ ! -f "$DATA" ]; then
    echo "ERROR: Data file not found: $DATA"
    exit 1
fi

echo "Merging:"
echo "  GEMC : $GEMC"
echo "  Data : $DATA"
echo "  -> $OUT"

hadd -f "$OUT" "$GEMC" "$DATA"

if [ $? -eq 0 ]; then
    echo "Done: $OUT"
else
    echo "ERROR: hadd failed"
    exit 1
fi
