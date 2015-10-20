#!/bin/bash

EXE=../../src/./mdff.x
sep1='======================================================================================================================'
sep2='---------------------------------------------------------'

echo $sep1
echo "EXAMPLE 4 : testing electric field gradient"
echo $sep1
echo ""
echo $sep2
echo "First example: Two charges in a box=10"
echo $sep2
echo " " 
cp TRAJFF.twocharges TRAJFF
cat > control.F << eof
&controltag
        calc='efg'
        lcoulomb  = .true.
        longrange = 'direct'
	cutlongrange = 1000.0D0
	lreduced=.true.
&end
&efgtag
        lefgprintall=.true.
        ncefg=1
&end
&fieldtag
	qch(1) =  1.0
	qch(2) = -1.0
        ncelldirect = 40
&end
&mdtag
&end
eof
echo $sep2
echo "Direct summation"
echo "ncelldirect = 40"
echo $sep2
echo ""
$EXE control.F > stdout_direct_twocharges
cat EFGALL > EFGALL.direct_twocharges
tail -n3 EFGALL 

echo $sep2
echo "Ewald summation"
echo $sep2
echo ""
cat > control.F << eof
&controltag
        calc='efg'
        lcoulomb  = .true.
        longrange = 'ewald'
	cutlongrange = 1000.0D0
	lreduced=.true.
&end
&efgtag
        lefgprintall=.true.
        ncefg=1
&end
&fieldtag
	qch(1) =  1.0
	qch(2) = -1.0
	kES = 11 11 11,
	alphaES = 0.8
&end
&mdtag
&end
eof
$EXE control.F > stdout_ewald_twocharges
cat EFGALL > EFGALL.ewald_twocharges
tail -n3 EFGALL 
echo $sep2
echo "From Nymand and Linse, J. Chem. Phys. 112, 6152 (2000)"
echo $sep2
echo ""
echo "Method           V1,xx                 V1,yy"
echo "ES             2.0003749            -1.0001874" 
echo "DS             2.0003749            -1.0001874"
echo ""
echo $sep1
echo $sep2
echo "Second example: cluster (9 atoms)"
echo "EFG Tensor principal components (NMR) :"
echo $sep2
echo ""
cp TRAJFF.cluster TRAJFF
echo $sep2
echo "Direct summation"
echo $sep2
echo ""
cat > control.F << eof
&controltag
        calc='efg'
        lcoulomb  = .true.
        longrange = 'direct'
	cutlongrange = 100.0d0
&end
&efgtag
        lefgprintall=.true.
        lefg_vasp_sign=.true.
	lefg_stat=.true.
        ncefg=1
	umin=-15.0
	smax=25.0
	vzzmin=-50.0
&end
&fieldtag
	qch(1) =  0.1 -0.8,
        ncelldirect = 40
&end
&mdtag
&end
eof
$EXE control.F > stdout_direct_cluster
cat NMRFF  > NMRFF.direct_cluster
tail -n10 NMRFF 
echo ""


echo $sep2
echo "Ewald summation"
echo $sep2
echo ""
cat > control.F << eof
&controltag
        calc='efg'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
&end
&efgtag
        lefgprintall=.true.
        lefg_vasp_sign=.true.
	lefg_stat=.true.
        ncefg=1
	umin=-15.0
	smax=25.0
	vzzmin=-50.0
&end
&fieldtag
	qch =  0.1 -0.8, 
	lautoES=.true.
	epsw=1e-8    
&end
&mdtag
&end
eof
$EXE control.F > stdout_ewald_cluster
cat NMRFF  > NMRFF.ewald_cluster
tail -n10 NMRFF 
echo ""
echo ""
echo $sep2
echo "GULP output "
echo $sep2
echo ""
echo "  -------------------------------------------------------------------------------"
echo "     Site no.    Diagonalised EFGs (V/Angs**2)          Asymmetry"
echo "                   xx        yy        zz               Parameter"
echo "  -------------------------------------------------------------------------------"
echo "         1       15.2454   16.4848  -31.7301             0.0391"
echo "         2       13.1974   15.5994  -28.7968             0.0834"
echo "         3       20.5185   23.5380  -44.0565             0.0685"
echo "         4       12.0418   14.6251  -26.6669             0.0969"
echo "         5        8.6231   10.8748  -19.4979             0.1155"
echo "         6       22.1626   24.1320  -46.2946             0.0425"
echo "         7        7.7041    9.9389  -17.6430             0.1267"
echo "         8       16.7300   17.0914  -33.8214             0.0107"
echo "         9       -0.8228   -3.4544    4.2772             0.6153"
echo "  -------------------------------------------------------------------------------"

#echo "                   EWALD(mdff)                      |                    DIRECT (mdff)                  |                       GULP" > COMP
#echo "       vxx          vyy          vzz          eta   |      vxx          vyy          vzz          eta   |      vxx          vyy          vzz          eta" >> COMP
#echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------">> COMP
#paste NMRFF.ewald_cluster NMRFF.direct_cluster gulp.out | tail -n9 | awk '{printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",$3,$4,$5,$7,$10,$11,$12,$14,$16,$17,$18,$20)}' >> COMP

#echo ""
#cat COMP
