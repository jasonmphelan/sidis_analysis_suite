#!/bash

./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots_no_match.root 1 0 &
./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots_no_match.root 0 0 &

#./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root 1 1 &
#./kaonAnalysis/makeKaonPlots /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root 0 1 &

wait

#./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_pi2k_AN.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_fits.root 2.5
#python ../plotting/plotKCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_fits.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/pi2k_fits
#python ../plotting/plotCorrections.py corrections_pi2k_AN.root pi2k

#./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_k2pi_AN.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_fits.root 2.5
#python ../plotting/plotKCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_fits.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/k2pi_fits
#python ../plotting/plotCorrections.py corrections_k2pi_AN.root k2pi


./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots_no_match.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_k2pi_AN_no_match.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_fits_no_match.root 2.5
#python ../plotting/plotCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/k2pi_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/k2pi_fits_no_match
python ../plotting/plotCorrections.py corrections_k2pi_AN_no_match.root k2pi

./kaonAnalysis/calcKaonCorr /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots_no_match.root /work/clas12/users/jphelan/sidis_analysis_suite/data/correctionFiles/corrections_pi2k_AN_no_match.root /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_fits_no_match.root 2.5
#python ../plotting/plotCutFits.py /work/clas12/users/jphelan/sidis_analysis_suite/histograms/analysis_note/pi2k_plots.root /volatile/clas12/users/jphelan/SIDIS/analysis_note/pi2k_fits_no_match
python ../plotting/plotCorrections.py corrections_pi2k_AN_no_match.root pi2k

./kaonAnalysis/kaonSys ../data/correctionFiles/corrections_pi2k_AN

python ../plotting/plotCorrections.py corrections_pi2k_AN_ratio.root pi2k
python ../plotting/plotCorrections.py corrections_k2pi_AN_ratio.root k2pi
