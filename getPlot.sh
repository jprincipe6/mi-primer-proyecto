sort -nk 2 score.fsc | awk '{print $2 "\t" $14}' > score_rmsd.dat
gnuplot rmsd_plot.p
