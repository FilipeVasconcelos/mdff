#!/bin/bash

EXEMDFF=mdff.x

cp config/POSFF.start .
rm data
for d in 1.30 1.25 1.20 1.175 1.15 1.125 1.10 1.08 1.06 1.04 1.02 1.0 0.95 0.9 0.8 0.7 0.6
do
    a3=`pc "256/$d"`
    a=`pc "$a3**(1./3.)"`
    mkdir -p runs/d$d
    cp POSFF.start runs/d$d/
    cp config/control_* runs/d$d/
    cd runs/d$d
    echo "256"         > POSFF
    echo "ex16 $d"    >> POSFF
    echo "$a 0.0 0.0" >> POSFF
    echo "0.0 $a 0.0" >> POSFF
    echo "0.0 0.0 $a" >> POSFF
    tail -n 260 POSFF.start >> POSFF
    $EXEMDFF control_md.F > log_md
    cp CONTFF ../../POSFF.start
    $EXEMDFF control_opt.F > log_opt
    cd ../../
    awk -v d=$d 'BEGIN{su=0.;sp=0.} {su+=$2;sp+=$4}END{print d,su/1000,sp/1000}' runs/d$d/ISTHFF >> data
done



