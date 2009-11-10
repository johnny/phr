set xlabel "m"
set ylabel "Zeit"
set logscale x

plot "stat13.dat" using 1:2 title 'stat13' with lines, \
   "stat102.dat" using 1:2 title 'stat102' with lines, \
   "dyn13.dat" using 1:2 title 'dyn13' with lines, \
   "dyn102.dat" using 1:2 title 'dyn102' with lines
