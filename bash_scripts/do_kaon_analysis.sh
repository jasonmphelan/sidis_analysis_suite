#!/bash

#./../build/kaonAnalysis/makeKaonPlots ../histograms/analysis_note/k2pi_plots_no_match.root 1 0 &
#./../build/kaonAnalysis/makeKaonPlots ../histograms/analysis_note/pi2k_plots_no_match.root 0 0 &

#./../build/kaonAnalysis/makeKaonPlots ../histograms/analysis_note/k2pi_plots.root 1 1 &
#./../build/kaonAnalysis/makeKaonPlots ../histograms/analysis_note/pi2k_plots.root 0 1 &

wait



./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/pi2k_plots.root ../data/correctionFiles/corrections_pi2k_AN.root ../histograms/analysis_note/pi2k_fits.root 2.5
#python ../plotting/plotKCutFits.py ../histograms/analysis_note/pi2k_plots.root ../histograms/analysis_note/pi2k_fits.root ../plotting/analysis_note/pi2k_fits
#python ../plotting/plotCorrections.py corrections_pi2k_AN.root pi2k

./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/k2pi_plots.root ../data/correctionFiles/corrections_k2pi_AN.root ../histograms/analysis_note/k2pi_fits.root 2.5
#python ../plotting/plotKCutFits.py ../histograms/analysis_note/k2pi_plots.root ../histograms/analysis_note/k2pi_fits.root ../plotting/analysis_note/k2pi_fits
#python ../plotting/plotCorrections.py corrections_k2pi_AN.root k2pi


./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/k2pi_plots_no_match.root ../data/correctionFiles/corrections_k2pi_AN_no_match.root ../histograms/analysis_note/k2pi_fits_no_match.root 2.5
#python ../plotting/plotCutFits.py ../histograms/analysis_note/k2pi_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/k2pi_fits_no_match
#python ../plotting/plotCorrections.py corrections_k2pi_AN_no_match.root k2pi

./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/pi2k_plots_no_match.root ../data/correctionFiles/corrections_pi2k_AN_no_match.root ../histograms/analysis_note/pi2k_fits_no_match.root 2.5
#python ../plotting/plotCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/pi2k_fits_no_match
#python ../plotting/plotCorrections.py corrections_pi2k_AN_no_match.root pi2k

./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/pi2k_plots.root ../data/correctionFiles/corrections_pi2k_AN_2sig.root ../histograms/analysis_note/pi2k_fits_2sig.root 2
./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/pi2k_plots.root ../data/correctionFiles/corrections_pi2k_AN_3sig.root ../histograms/analysis_note/pi2k_fits_3sig.root 3

./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/k2pi_plots.root ../data/correctionFiles/corrections_k2pi_AN_2sig.root ../histograms/analysis_note/k2pi_fits.root 2
./../build/kaonAnalysis/calcKaonCorr ../histograms/analysis_note/k2pi_plots.root ../data/correctionFiles/corrections_k2pi_AN_3sig.root ../histograms/analysis_note/k2pi_fits.root 3



./../build/kaonAnalysis/kaonSys ../data/correctionFiles/corrections_pi2k_AN
./../build/kaonAnalysis/kaonSys ../data/correctionFiles/corrections_k2pi_AN

python ../plotting/plotCorrections.py corrections_pi2k_AN_ratio.root pi2k
python ../plotting/plotCorrections.py corrections_k2pi_AN_ratio.root k2pi
