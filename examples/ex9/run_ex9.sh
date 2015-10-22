#!/bin/bash

pwds=$PWD
EXEMDFF=mdff.x
EXECP2K=cp2k

echo "generate control files"
cd config 
./gen_tests
ntest=`wc -l DOALLC | awk '{print $1}'`
cd ..

for ((k=1;k<=$ntest;k++))
do
	echo " ===================================================================== "
	rep=test.$k
	echo " doing $rep ..."
	# ----------------------------------------
	echo " cp2k running"
	if [ -d cp2k/$rep ]; then
		echo " cp2k/$rep already exists... cleaning"
		rm cp2k/$rep/*
	else
		mkdir -p cp2k/$rep
	fi
	cd cp2k/$rep
	cp $pwds/config/run_cp2k run
	cp $pwds/config/pot_cp2k/*CP2K.POT .
	cp $pwds/config/pot_cp2k/input_$rep.cp2k input.cp2k
	./run $EXECP2K
	cd $pwds
	# ----------------------------------------
	echo " mdff running"
	if [ -d mdff/$rep ]; then
		echo " mdff/$rep already exists... cleaning"
		rm mdff/$rep/*
	else
		mkdir -p mdff/$rep
	fi
	cd mdff/$rep
	cp $pwds/config/POSFF .
	cp $pwds/config/run_mdff run
	cp $pwds/config/pot_mdff/control_$rep.F control.F
	./run $EXEMDFF $rep
	cd $pwds
	# ----------------------------------------

done
