#!/bin/bash
# files :
#  		  config/control.F  (MDFF)
#		  config/POSFF      (MDFF)
#		  config/CONFIG     (DLPOLY)
#		  config/CONTROL    (DLPOLY)
#		  config/FIELD      (DLPOLY)
#                 config/input.cp2k (CP2K)
#
# ====================================================
# user settings

EXECP2K=cp2k
EXEDLPOLY=DLPOLY.X
EXEMDFF=mdff.x

do_cp2k=true
do_dlpoly=true

# =====================================================

sep="=================================================="
echo $sep
echo "# Example 0: Minimum settings with FCC structure"
echo "# FCC LJ a=1.5414"
echo "# natm = 256"
echo "# more details in control.F file"

rm -rf mdff/
mkdir mdff
cd mdff
echo 
echo "# mdff calculation "
cp ../config/control.F .
cp ../config/POSFF .
$EXEMDFF control.F > stdout 
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
	cd ..
fi

if $do_cp2k; then
	echo "# cp2k calculation requested "
	rm -rf cp2k 
	mkdir cp2k
	cd cp2k
	cp ../config/input.cp2k .
	$EXECP2K -i input.cp2k -o output.cp2k
	cd ..
fi


echo
echo $sep
echo 
# print energy informations
etot=`grep Etot mdff/stdout | awk '{print $NF}'` 
etot2=`echo ${etot} | sed -e 's/[eE]+*/\\*10\\^/'`
natm=`grep natm mdff/stdout | awk '{print $NF}'`
lr=`grep "long range correction (init) energy" mdff/stdout |  awk '{print $NF}'`
etot3=`echo "scale=8;$etot2" | bc -l | awk '{printf("%.4f\n", $1)}'`
etotlr=`echo "scale=8;($etot2)+($lr)" | bc -l | awk '{printf("%.4f\n", $1)}'`
echo "Etot (mdff)                  = " $etot3
echo "Etot + long range correction = " $etotlr
echo ""


if [ $do_cp2k = "true" ] && [ -f cp2k/FCC_LJ-1.ener ]
then
	etotcp2k=`tail -n 1 cp2k/FCC_LJ-1.ener | awk '{printf ("%.4f\n", $5*27.211385)}'`
	echo "cp2k                         = " $etotcp2k 
fi
if [ $do_dlpoly = "true" ] && [ -f dl_poly/STATIS ]
then
	etotdl=`head -n 4 dl_poly/STATIS | tail -n 1 | awk '{printf("%.4f\n", $1*1.0)}'`
	echo "dl_poly                      = " $etotdl 
fi


exit 0
