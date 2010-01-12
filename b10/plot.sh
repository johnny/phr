nop=(2 4 6 8 10 12)
i=3
plot=""
for procs in ${nop[@]}
do
    plot="$plot \"ue26.dat\" u 1:(\$$i/\$2) title \"$procs threads\" with linespoints,"
    i=$(($i + 1 ))
done
plot=${plot%\,}
echo "set xlabel \"Bodies\"
set ylabel \"Speedup\"
set term png
set out \"plot.png\"
plot $plot" | gnuplot -persist

exit 0

