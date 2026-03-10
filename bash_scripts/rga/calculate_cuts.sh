#!/bash

./../../build/cutAnalysis/fiducials ../../histograms/rga/fiducials.root 10.6 1 &
./../../build/cutAnalysis/calcSFCuts ../../histograms/rga/SF_cuts.root 10.6 1 &
./../../build/cutAnalysis/calcVertex ../../histograms/rga/vertex_cuts.root 10.6 1 &