#!/bin/bash

echo "#!/usr/bin/gnuplot " > plot.dos
echo "set term postscript enhanced eps color 14" >> plot.dos
echo 'set output "dos.eps"' >> plot.dos

#echo "set style line 1 lw 2 lt 1 lc 1" >> plot.dos
#echo "set style line 2 lw 2 lt 1 lc 2" >> plot.dos
#echo "set style line 3 lw 2 lt 1 lc 3" >> plot.dos
echo "set style line 1 lw 2 lt 1 " >> plot.dos
echo "set style line 2 lw 2 lt 2 " >> plot.dos
echo "set style line 3 lw 2 lt 3 " >> plot.dos
echo "set style line 4 lw 2 lt 4 " >> plot.dos
echo "set style line 5 lw 2 lt 5 " >> plot.dos

echo "set xlabel '{/Symbol w}'" >> plot.dos
echo "set ylabel 'g({/Symbol w})'" >> plot.dos

echo "p 'DOSKFF.vib+dos' u 1:2 w l title 'calc=vib+dos' ls 1 ,'' u 1:3 w l title 'T1' ls 2,'' u 1:4 w l title 'T2' ls 3,'' u 1:5 w l title 'L' ls 4,'DOSFF.vib' w l title 'calc=vib' ls 5">> plot.dos
echo "  #    EOF" >> plot.dos
chmod u+x plot.dos
./plot.dos
gv dos.eps

