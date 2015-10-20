#!/bin/bash

echo "#!/usr/bin/gnuplot -persist" > plot.band
echo "set term postscript enhanced eps color 20" >> plot.band
echo 'set output "band.eps"' >> plot.band

echo "set style line 1 lw 2 lt 1 lc 1" >> plot.band
echo "set style line 2 lw 2 lt 1 lc 2" >> plot.band
echo "set style line 3 lw 2 lt 1 lc 3" >> plot.band

echo "set ylabel 'energy {/Symbol w}'" >> plot.band
echo 'set format x ""' >> plot.band
echo "set xtics ('{/Symbol G}' 0,'{/Symbol G}' -2,'X' -1,'L' 1)" >> plot.band

echo "p 'DOSKFF.vib+band_LtoG' u 2:5 w l title 'T1' ls 1 ,\
        'DOSKFF.vib+band_LtoG' u 2:6 w l title 'T2' ls 2 ,\
        'DOSKFF.vib+band_LtoG' u 2:7 w l title 'L'  ls 3 ,\
        'DOSKFF.vib+band_GtoX' u (-\$4):5 w l title '' ls 1 ,\
        'DOSKFF.vib+band_GtoX' u (-\$4):6 w l title '' ls 2 ,\
        'DOSKFF.vib+band_GtoX' u (-\$4):7 w l title '' ls 3 ,\
        'DOSKFF.vib+band_XtoG' u (\$4-2):5 w l title '' ls 1 ,\
        'DOSKFF.vib+band_XtoG' u (\$4-2):6 w l title '' ls 2 ,\
        'DOSKFF.vib+band_XtoG' u (\$4-2):7 w l title '' ls 3  ">> plot.band
echo "  #    EOF" >> plot.band
chmod u+x plot.band
./plot.band
gv band.eps

