#!/bin/bash
# files :
#                 config/control.F  (MDFF)
#                 config/POSFF      (MDFF)
#                 config/CONFIG     (DLPOLY)
#                 config/CONTROL    (DLPOLY)
#                 config/FIELD      (DLPOLY)
#                 config/input.cp2k (CP2K)
#
# ====================================================
# user settings

EXEDLPOLY=DLPOLY.X
EXEMDFF=mdff.x
bin=/home/filipe/dev/mdff/bin

do_dlpoly=true

# =====================================================

sep="=================================================="
echo $sep
echo "#Example 1: LJ (argon) fcc structure at low temperature"
echo "The configuration is readed in POSFF"
echo "more info in control.F and stdout"

rm -rf mdff/
mkdir mdff
cd mdff
echo 
echo "# mdff calculation "
cp ../config/control.F .
cp ../config/POSFF .
$EXEMDFF control.F > stdout
${bin}/poszi.py -i OSZIFF -n
cd ..

if $do_dlpoly; then
        echo "# dlpoly calculation requested "
        rm -rf dl_poly
        mkdir dl_poly
        cd dl_poly
        cp ../config/CONTROL .
        cp ../config/FIELD .
        cp ../config/CONFIG .
        $EXEDLPOLY
	${bin}/read_statis > ene
        cd ..
fi

lr=`grep "long range correction (init) energy" mdff/stdout |  awk '{print $NF}'`
echo $lr
cat > plot << eof
#!/usr/bin/gnuplot 
reset
set term x11
set title "Total energy MD ex1."
set xlabel "time"
set ylabel "E_{tot}"
p 'dl_poly/ene' u 1:2 w l title "DL_POLY"
pause -1
eof
chmod u+x plot
./plot
