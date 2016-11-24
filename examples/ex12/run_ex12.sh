#!/bin/bash

EXEDLPOLY=cp2k
EXEMDFF=mdff.x

cp2k=true
mdff=false

if $cp2k; then
	echo "running CP2K test (ex12)"
	rm -rf cp2k
	mkdir -p cp2k
	cd cp2k
	cp ../config/CONFIG .
	cp ../config/CONTROL .
	cp ../config/FIELD .
		$EXEDLPOLY 
		tail -n +4 RDFDAT > RDFDAT.$ensemble
	cd ..
fi

if $mdff; then
	echo "running mdff test (ex8)"
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
p 'cp2k/RDFDAT.$ensemble' w l ,'mdff/GRTFF.$ensemble' w l
eof
chmod u+x plot.gr
./plot.gr

