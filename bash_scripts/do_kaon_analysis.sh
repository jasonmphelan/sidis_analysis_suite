#!/bash

./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root 1 &
./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root 0 &

wait

./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_k2pi_AN.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_fits.root
python ./plotting/plotCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/k2pi_fits
python ../plotting/plotCorrections.py corrections_k2pi_AN.root k2pi

./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_pi2k_AN.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_fits.root
python ./plotting/plotCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/pi2k_fits
python ../plotting/plotCorrections.py corrections_pi2k_AN.root pi2k

