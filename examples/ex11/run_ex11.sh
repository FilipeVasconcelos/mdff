#!/bin/bash

pwds=$PWD


dir=random_structs

rm -rf $dir
mkdir -p $dir
cp config/control_stochio.F $dir/
cd $dir

random_struct+stochio -s 102 -f 129

cd $pwds 


