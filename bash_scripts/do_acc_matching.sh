#!/bash

#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_1.root 0 .01 &

#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_5.root 0 .05 &
./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_map.root 0 .1 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_15.root 0 .15 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_20.root 0 .20 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_25.root 0 .25 &
#./../build/cutAnalysis/acceptanceMatching ../data/acceptance_matching/matchCut2D_50.root 0 .50 &

wait

python ../plotting/plotAcceptanceMatching.py ../data/acceptance_matching/matchCut2D_map.root ../plotting/analysis_note/acceptance_matching
python ../plotting/plotAcceptanceMatching.py ../data/acceptance_matching/matchCutPi2K_map.root ../plotting/analysis_note/acceptance_matching
python ../plotting/plotAcceptanceMatching.py ../data/acceptance_matching/matchCutK2Pi_map.root ../plotting/analysis_note/acceptance_matching
