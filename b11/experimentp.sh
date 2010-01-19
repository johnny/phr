#!/bin/bash

nop=(2 4 8 16 32)
iter=100
rm ue28.dat

for N in $*
do
    echo Number Of Bodies: $N
    out="$N"

    for procs in ${nop[@]}
    do
        echo running synchronous calculation on $procs procs
        mpirun -machinefile mpihosts -np $procs nbody_mpi $N $iter 1 > ps.$N.$procs
        sum=0
        while read num
        do
            sum=`echo $sum + $num | bc`
        done < <(cat ps.$N.$procs)
        scaled=`echo "scale=4 ; $sum / $iter" | bc`
        out="$out $scaled"

        echo running asynchronous calculation on $procs procs
        mpirun -machinefile mpihosts -np $procs nbody_mpi_a $N $iter 1 > pa.$N.$procs
        sum=0
        while read num
        do
            sum=`echo $sum + $num | bc`
        done < <(cat pa.$N.$procs)
        scaled=`echo "scale=4 ; $sum / $iter" | bc`
        out="$out $scaled"
    done

    echo $out >> ue28.dat
done

i=2
plot=""
for procs in ${nop[@]}
do
    plot="$plot \"ue28.dat\" u 1:(\$$(($i + 1 ))/\$$i) title \"$procs threads\" with linespoints,"
    i=$(($i + 2 ))
done
plot=${plot%\,}
echo "set xlabel \"Bodies\"
set ylabel \"Speedup\"
set term png
set out \"plot.png\"
plot $plot" | gnuplot -persist

exit 0
