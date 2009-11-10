set xlabel "N"
set ylabel "MFlops"

plot "plot.data" using 1:2 title 'gauss' with lines, \
   "plot.data" using 1:3 title 'gauss blocked'with lines
