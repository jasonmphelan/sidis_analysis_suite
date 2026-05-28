#!/bin/bash
# Run tight_reskim over final_skims/10.X (data) and final_skims/GEMC (MC)
# Usage: ./bash_scripts/run_tight_reskim.sh  (from project root)

BINARY=./../build/skimmers/tight_reskim
IN_DATA=./../trees/final_skims
IN_MC=./../trees/final_skims/GEMC
OUT_DATA=./../trees/final_skims
OUT_MC=./../trees/final_skims/GEMC
TARGET=0 # 1 = RGA/proton

# Data (RunType=0)
#$BINARY $IN_DATA/10.2/final_skim.root $OUT_DATA/10.2/tight_skim_new.root 10.2 $TARGET 0 &
#$BINARY $IN_DATA/10.4/final_skim.root $OUT_DATA/10.4/tight_skim_new.root 10.4 $TARGET 0 &
#$BINARY $IN_DATA/10.6/final_skim.root $OUT_DATA/10.6/tight_skim_new.root 10.6 $TARGET 0 &

# MC (RunType=1)
$BINARY $IN_MC/final_skim_10.2_merged.root $OUT_MC/tight_skim_10.2.root 10.2 $TARGET 1 &
$BINARY $IN_MC/final_skim_10.4.root $OUT_MC/tight_skim_10.4.root 10.4 $TARGET 1 &
$BINARY $IN_MC/final_skim_10.6.root $OUT_MC/tight_skim_10.6.root 10.6 $TARGET 1 &

wait
echo "Done."
