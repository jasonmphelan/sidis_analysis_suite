#!/bash

bash do_kaons.sh
wait
bash make_ratios.sh
wait
bash make_ratio_plots.sh
