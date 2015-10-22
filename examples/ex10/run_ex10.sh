#!/bin/bash
# ====================================================

sep="======================================================"
echo $sep
echo "# Example 10 : Generating control.F files for a set"
echo "               of runs. see README for more details"
echo $sep

dir=control_files
rm -r $dir
mkdir -p $dir
cp config/* $dir
cd $dir

gen_quench_controls

