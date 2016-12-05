#!/bin/bash

EXEMDFF=mdff.x

echo "EXAMPLE 3 : "

rep=runs
if [ -d $rep ] 
then
	rm -r $rep 
else
	mkdir -p $rep
fi
cd $rep
	cp ../config/control_opt.F .
	cp ../config/run.py .
	./run.py
cd ..

