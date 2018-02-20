#!/bin/bash

EXEMDFF=mdff.x

cp config/POSFF.start .
rm data
for d in 1.30 1.25 1.20 1.175 1.15 1.125 1.10 1.09 1.08 1.07 1.06 1.05 1.04 1.02 1.0 0.975 0.95 0.9 0.85 0.8 0.7 0.6
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
    mpirun -np 2 $EXEMDFF control_md.F > log_md
    cp CONTFF ../../POSFF.start
    mpirun -np 2 $EXEMDFF control_opt.F > log_opt
    cd ../../
    awk -v d=$d 'BEGIN{su=0.;sp=0.} {su+=$2;sp+=$4}END{print d,su/1000,sp/1000}' runs/d$d/ISTHFF >> data
done

cat > plot << eof
#!/usr/bin/gnuplot
reset
set terminal x11 size 300,600
set termoption enhanced
set encoding utf8
set xr[0.6:1.3]
set pointsize 3.0
set size 1.0,1.0
set multiplot layout 2, 1 title "Figure 1 of PRB 68, 144202 (2003)";
set size 1.0,0.45
set xlabel "{/symbol r}"
set ylabel "<U>"
p 'data' u 1:2 w lp pt 12 title ""
set size 1.0,0.45
set xlabel "{/symbol r}"
set ylabel "<P>"
p 'data' u 1:3 w lp pt 12 title ""
unset multiplot
pause -1
eof
chmod u+x plot
./plot


