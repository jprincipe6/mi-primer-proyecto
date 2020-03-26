
set xlabel "rmsd"
set ylabel "score"
set term png
set output "score_rmsd_plot.png"
plot "score_rmsd.dat" u 2:1 title "rmsd/score"
replot
#pause -1 "continuar"
