#!/bin/bash

rm mpihosts

n=0
for i in `seq -w 50`; do
    load=`ssh "cip$i" cut -f 1 -d \" \" /proc/loadavg 2>/dev/null|tr -d .`
    if test -n "$load" && test "$load" -lt 50; then
        echo "Adding cip$i..."
        echo "cip$i" >> mpihosts
        n=$(($n+1))
    fi
done

echo "$n hosts added."
