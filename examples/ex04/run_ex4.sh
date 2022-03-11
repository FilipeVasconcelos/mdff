#!/bin/bash

EXE=mdff.x
sep1='======================================================================================================================'
sep2='---------------------------------------------------------'

echo $sep1
echo "EXAMPLE 4 : testing electric field gradient / NMR"
echo $sep1
echo ""
echo $sep2
echo "First example: From [1] Nymand and Linse, J. Chem. Phys. 112, 6152 (2000) "
echo "               Two charges / two dipoles / Two charges and two polarizabilities (reduced units)"
echo $sep2
echo " " 
cat > control_twocharges.F << eof
! field calculation
&controltag
        calc='md'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
	lreduced=.true.            ! reduced units 
	lreducedN=.false.          ! NOT REDUCED BY N
&end
&fieldtag
	qch =  1.0 -1.0 , 
	kES = 11 11 11,            ! Parameters from [1]
	alphaES = 0.8              ! from [1]

        doefg=.true.               ! do EFG
        doefield=.true.            ! do Electric Field
        lwrite_efg=.true.          ! and write ...
        lwrite_ef=.true.
&end
&mdtag
&end
eof
cat > control_twodipoles.F << eof
! field calculation
&controltag
        calc='md'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
	lreduced=.true.            ! reduced units 
	lreducedN=.false.          ! NOT REDUCED BY N
&end
&fieldtag
        dip(1,1) =  1.0 
        dip(2,1) =  1.0 
	kES = 11 11 11,            ! Parameters from [1]
	alphaES = 0.8              ! from [1]

        doefg=.true.               ! do EFG
        doefield=.true.            ! do Electric Field
        lwrite_efg=.true.          ! and write ...
        lwrite_ef=.true.
&end
&mdtag
&end
eof

cat > control_twocharges_twopola.F << eof
! field calculation
&controltag
        calc='md'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
	lreduced=.true.            ! reduced units 
	lreducedN=.false.          ! NOT REDUCED BY N
	lsurf=.true.
&end
&fieldtag
	qch =  1.0 -1.0 , 
	kES = 11 11 11,            ! Parameters from [1]
	alphaES = 0.8              ! from [1]
        min_scf_pol_iter = 3
        max_scf_pol_iter = 3
        conv_tol_ind = 1e-6
	!algo_moment_from_pola = 'scf'
	algo_moment_from_pola = 'scf_kO_v2'


	lpolar(3) =.true.
	pol(3,1,1) = 0.1d0

	lpolar(4) =.true.
	pol(4,1,1) = 0.1d0
	pol(4,2,2) = 0.1d0
	pol(4,3,3) = 0.1d0
	ldip_damping = .false.

        doefg=.true.               ! do EFG
        doefield=.true.            ! do Electric Field
        lwrite_efg=.true.          ! and write ...
        lwrite_ef=.true.
&end
&mdtag
&end
eof


tag=twocharges
cp POSFF.$tag POSFF 
$EXE control_$tag.F > stdout_$tag
cat EFGALL > EFGALL.$tag
cat EFALL > EFALL.$tag
cat CONTFF > CONTFF.$tag


Utmp=`grep Utot stdout_$tag | awk '{print $NF}'`
Utmp2=`echo ${Utmp} | sed -e 's/[eE]+*/\\*10\\^/'`
U=`echo "scale=8;$Utmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Etmp=`tail -n+10 EFALL | head -n 1 | awk '{print $3}'`
Etmp2=`echo ${Etmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Ex=`echo "scale=8;$Etmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $3}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vxx=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $4}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vyy=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Ftmp=`tail -n+10 CONTFF | head -n 1 | awk '{print $8}'`
Ftmp2=`echo ${Ftmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Fx=`echo "scale=8;$Ftmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`


echo $sep2
echo "From Nymand and Linse, J. Chem. Phys. 112, 6152 (2000)"
echo $sep2
echo "                                        Two Charges "
echo "Method           U              E1,x              V1,xx          V1,yy            f1,x"
echo "ES             -1.0021255      0.9956865        2.0003749       -1.0001874       0.9956865"   
echo ""
echo "MDFF(this run) $U     $Ex       $Vxx      $Vyy      $Fx"
echo ""

tag=twodipoles
cp POSFF.$tag POSFF 
$EXE control_$tag.F > stdout_$tag
cat EFGALL > EFGALL.$tag
cat EFALL > EFALL.$tag
cat CONTFF > CONTFF.$tag

Utmp=`grep Utot stdout_twodipoles | awk '{print $NF}'`
Utmp2=`echo ${Utmp} | sed -e 's/[eE]+*/\\*10\\^/'`
U=`echo "scale=8;$Utmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Etmp=`tail -n+10 EFALL | head -n 1 | awk '{print $3}'`
Etmp2=`echo ${Etmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Ex=`echo "scale=8;$Etmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $3}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vxx=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $4}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vyy=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Ftmp=`tail -n+10 CONTFF | head -n 1 | awk '{print $8}'`
Ftmp2=`echo ${Ftmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Fx=`echo "scale=8;$Ftmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`


echo "                                        Two Dipoles"
echo "ES             -2.0087525      2.0087525        5.9992460       -2.9996230       5.9992460 "   
echo ""
echo "MDFF(this run) $U     $Ex       $Vxx      $Vyy      $Fx"


tag=twocharges_twopola
cp POSFF.$tag POSFF 
$EXE control_$tag.F > stdout_$tag
cat EFGALL > EFGALL.$tag
cat EFALL > EFALL.$tag
cat CONTFF > CONTFF.$tag


Utmp=`grep Utot stdout_$tag | awk '{print $NF}'`
Utmp2=`echo ${Utmp} | sed -e 's/[eE]+*/\\*10\\^/'`
U=`echo "scale=8;$Utmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Etmp=`tail -n+10 EFALL | head -n 1 | awk '{print $3}'`
Etmp2=`echo ${Etmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Ex=`echo "scale=8;$Etmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $3}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vxx=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Vtmp=`tail -n+10 EFGALL | head -n 1 | awk '{print $4}'`
Vtmp2=`echo ${Vtmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Vyy=`echo "scale=8;$Vtmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`

Ftmp=`tail -n+10 CONTFF | head -n 1 | awk '{print $8}'`
Ftmp2=`echo ${Ftmp} | sed -e 's/[eE]+*/\\*10\\^/'`
Fx=`echo "scale=8;$Ftmp2" | bc -l | awk '{printf("%.8f\n", $1)}'`


echo "                                        Two charges and two polarizabilities"
echo "ES             -1.012 "   
echo ""
echo "MDFF(this run) $U     $Ex       $Vxx      $Vyy      $Fx"




echo $sep1
echo $sep2
echo "Second example: cluster (9 atoms)"
echo "EFG Tensor principal components (NMR) :"
echo $sep2
echo ""



cp POSFF.cluster POSFF 


echo ""
cat > control_1.F << eof
&controltag
        calc='md'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
&end
&fieldtag
	qch =  0.1 -0.8, 
	lautoES=.true.
	epsw=1e-12
        doefg=.true.    
        doefield=.true.    
        lwrite_efg=.true.
        lwrite_ef=.true.
&end
&mdtag
&end
eof
cat > control_2.F << eof
&controltag
        calc='efg'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
&end
&efgtag
        lefg_restart=.true.
        lefg_stat=.true.
        lefg_vasp_sign=.true.
        ncefg=1
	umin=-15.0
	smax=25.0
	vzzmin=-50.0
&end
&fieldtag
&end
&mdtag
&end
eof
$EXE control_1.F > stdout_cluster_1
$EXE control_2.F > stdout_cluster_2
cat NMRFF  > NMRFF.cluster
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
