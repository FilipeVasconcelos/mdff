#!/bin/bash

EXEMDFF=mdff.x

echo "EXAMPLE 3 : "
echo "Global optimisation of clusters"
echo "Random structures are generated"
echo "and Local minima ( m1qn3 or lbfgs ) are found for each of them"
echo "better method as to be used for larger cluster"

rep=runs
if [ -d $rep ] 
then
	rm -r $rep 
	mkdir -p $rep
else
	mkdir -p $rep
fi
cd $rep
	cp ../config/control_opt.F .
	cp ../config/run.py .
	./run.py
cd ..

