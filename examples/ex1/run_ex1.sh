#!/bin/bash
echo "#Example 1: LJ (argon) fcc structure at low temperature"
echo "The configuration is readed in POSFF"
echo "more info in control.F and stdout"
EXE=../../bin/mdff.x 
mpirun -n 2 $EXE control.F > stdout

lr=`grep "long range correction (init) energy" stdout |  awk '{print $NF}'`
echo $lr
poszi -i OSZIFF -n
cat > plot << eof
#!/usr/bin/gnuplot -persist
reset
set term x11
set title "Total energy MD ex2."
set xlabel "time"
set ylabel "E_{tot}"
p 'dl_poly/ene' u 1:2 w l title "DL_POLY",'etot_l' u 6:(\$9+($lr)) w l title "MDFF"
eof
chmod u+x plot
./plot
