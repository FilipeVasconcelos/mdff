#!/bin/bash

EXEDLPOLY=DLPOLY.X
EXEMDFF=mdff.x

#dlpoly=true
mdff=true
ensemble=nve

#if $dlpoly; then
#	echo "running DL_POLY test (ex7)"
#	rm -rf dlpoly
#	mkdir -p dlpoly
#	cd dlpoly
#		cp ../config/CONFIG .
#		cp ../config/FIELD .
#		cp ../config/CONTROL .
#		$EXEDLPOLY
#		tail -n +4 RDFDAT > RDFDAT.$ensemble
#	cd ..
#fi

if $mdff; then
	echo "running mdff test (ex7)"
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
#!/usr/bin/gnuplot -persist
reset
p 'dlpoly/RDFDAT.$ensemble' w l ,'mdff/GRTFF.$ensemble' w l
eof
chmod u+x plot.gr
./plot.gr

