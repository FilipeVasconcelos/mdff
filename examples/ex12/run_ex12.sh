#!/bin/bash

EXECP2K=cp2k
EXEMDFF=mdff.x

cp2k=true
mdff=true

if $cp2k; then
	echo "running CP2K test (ex12)"
	rm -rf cp2k
	mkdir -p cp2k
	cd cp2k
	cp ../config/mdff/control_gr.F .
	cp ../config/cp2k/input_pim.cp2k .
	cp ../config/cp2k/*.POT .
	$EXECP2K -i input_pim.cp2k -o output_pim.cp2k
	mv SiO2-1.cell SiO2-pos-1.cell
	xyztoff -i SiO2-pos-1.xyz -o TRAJFF
	mpirun -np 2 $EXEMDFF control_gr.F > stdout_gr

	cd ..
fi

if $mdff; then
	echo "running mdff test (ex12)"
	rm -rf mdff 
	mkdir -p mdff
	cd mdff
	cp ../config/mdff/control_*.F .
	cp ../config/mdff/POSFF .
	mpirun -np 2 $EXEMDFF control_md.F > stdout_md
	$EXEMDFF control_gr.F > stdout_gr
	cd ..
fi

cat > plot.gr << eof
#!/usr/bin/gnuplot 
reset
p 'cp2k/GRTFF' w l ,'mdff/GRTFF' w l
pause -1
eof
chmod u+x plot.gr
./plot.gr

