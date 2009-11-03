set xlabel "Blockgroesse"
set ylabel "MFlops"

plot "plotc.data" using 1:2 title 'gauss blocked, Problemgroesse 128' with lines
