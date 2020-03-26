set xlabel "nstructs"
set ylabel "Time (ms)"
set term png
set output "time_plot.png"
set style line 1
set linetype 1 linewidth 2
plot "salida_time.txt" u ($0+1):1 with linespoints ls 1 title "Time for Nstructs"
replot
#pause -1 "continuar"
