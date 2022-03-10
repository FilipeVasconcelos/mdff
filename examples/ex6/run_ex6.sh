#!/bin/bash

EXEDLPOLY=DLPOLY.X
EXEMDFF=mdff.x

dlpoly=true
mdff=true
ensemble=nvt

if $dlpoly; then
	echo "running DL_POLY test (ex6)"
	rm -rf dlpoly
	mkdir -p dlpoly
	cd dlpoly
	cp ../config/CONFIG .
	cp ../config/CONTROL .
	cp ../config/FIELD .
		$EXEDLPOLY
		tail -n +4 RDFDAT > RDFDAT.$ensemble
	cd ..
fi


if $mdff; then
	echo "running mdff test (ex6)"
	rm -rf mdff
       	mkdir -p mdff	
	cd mdff
	cp ../config/control_* .
	cp ../config/POSFF .
	$EXEMDFF control_md.F > stdout_md.$ensemble
	$EXEMDFF control_gr.F > stdout_gr.$ensemble
	mv GRTFF GRTFF.$ensemble
	cd ..
fi

cat > plot.gr << eof
#!/usr/bin/gnuplot 
reset
p 'dlpoly/RDFDAT.$ensemble' w l ,'mdff/GRTFF.$ensemble' w l
pause -1
eof
chmod u+x plot.gr
./plot.gr

