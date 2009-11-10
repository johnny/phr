set xlabel "m"
set ylabel "Zeit"
set logscale x

plot "plota128.data" using 1:3 title '128' with lines, \
   "plota256.data" using 1:3 title '256'with lines,\
   "plota512.data" using 1:3 title '512'with lines
