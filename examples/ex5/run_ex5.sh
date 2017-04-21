#!/bin/bash

EXE=/home/filipe/dev/mdff/bin/mdff.x
echo "===================================================="
echo "EXAMPLE 5 : Density of states of FCC Lennard-Jones"
echo "===================================================="
echo ""
echo "----------------------------------"
echo "Hessian matrix and dos calculation"
echo "----------------------------------"
$EXE control_vib.F > stdout.vib
mv DOSFF DOSFF.vib
#./plot_vectff.sh 15.0 0.01 256
echo ""
echo "---------------------------------------------"
echo "complete DOS "
echo "(nkphon+1) x (nkphon+1) x (nkphon+1) k-points"
echo "generated k-points in IBZKPTFF" 
echo "---------------------------------------------"
mpirun -np 2 $EXE control_vib+dos.F > stdout.vib+dos
mv DKFF DKFF.vib+dos
mv DOSKFF DOSKFF.vib+dos
echo "plot dos"
./plot_dos_ex5.sh 
echo ""
echo "---------------------------------------------"
echo "Band calculation"
echo "---------------------------------------------"
echo "from L to Gamma"
$EXE control_vib+band_LtoG.F > stdout.vib+band_LtoG
mv DOSKFF DOSKFF.vib+band_LtoG
echo "from Gamma to X"
$EXE control_vib+band_GtoX.F > stdout.vib+band_GtoX
mv DOSKFF DOSKFF.vib+band_GtoX
echo "from X to Gamma"
$EXE control_vib+band_XtoG.F > stdout.vib+band_XtoG
mv DOSKFF DOSKFF.vib+band_XtoG
echo ""
echo "plot band stucture"
./plot_band_ex5.sh 
echo "end of example 5"
