#!/bin/bash

nop=(2 4 8 16 32)
iter=100
rm ue26.dat

for N in $*
do
    echo Number Of Bodies: $N
    out="$N"

    echo running vanilla calculation
    ./nbody_vanilla $N $iter 1 > s.$N
    sum=0
    while read num
    do
        sum=`echo $sum + $num | bc`
    done < <(cat s.$N)
    scaled=`echo "scale=4 ; $sum / $iter" | bc`
    out="$out $scaled"
    
    for procs in ${nop[@]}
    do
        echo running parallell calculation on $procs procs
        mpirun -machinefile mpihosts -np $procs nbody_mpi $N $iter 1 > p.$N.$procs
        sum=0
        while read num
        do
            sum=`echo $sum + $num | bc`
        done < <(cat p.$N.$procs)
        scaled=`echo "scale=4 ; $sum / $iter" | bc`
        out="$out $scaled"
    done

    echo $out >> ue26.dat
done
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
