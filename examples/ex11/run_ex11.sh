#!/bin/bash

pwds=$PWD


dir1=random_structs 
dir2=md_runs

sep="======================================================"
echo $sep
echo "# Example 11 : Generating 5 random configuration of SiO2 with"
echo "               102  atoms and execute  a short  equilibration"
echo "               run for each of them."
echo "               ( NVE ~ 2 ps Rigid Ions and Polarizable Ions models ) "
echo $sep

rm -rf $dir1
mkdir -p $dir1
cp config/control_stochio.F $dir1/
cd $dir1

for i in {1..5}
do
	random_struct+stochio -s 102 -f 102 -t $i
done

cd $pwds 


rm -rf $dir2
mkdir -p $dir2
cp config/*.POT $dir2
cp config/quench.config $dir2
cp config/control_template.F $dir2
cp $dir1/POSFF.* $dir2
cd $dir2

gen_quench_controls

for file in POSFF.randomPY.*
do
	rm -rf r$file
	mkdir -p r$file
	cd r$file
	cp ../$file POSFF.start
	cp ../control_*.F .
	cp ../DOALLMD .
	cp ../quench.config .

	#DoAll_mdff.sh --para=2
	DoAll_mdff.sh --serial

	cd ..

done

