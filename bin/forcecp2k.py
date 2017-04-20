#!/usr/bin/env python
# ===============================================================
# date       : 1er decembre 2016 
# author     : fmv
# decription : extract forces from cp2k output 
# ===============================================================

import datetime
import argparse
from constants import float_format
import re                                 

def main_parser():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i','--input', help='Input file name (CP2K output)',required=True,dest="input_file")
    parser.add_argument('-o','--output',help='Output file name (forces)',required=True,dest="output_file")

    args = parser.parse_args()

    return args




if __name__ == "__main__":

    separator=60*"="
    now = datetime.datetime.now()

    print separator
    print now.strftime("%Y-%m-%d %H:%M")
    print "author : filipe.manuel.vasconcelos@gmail.com"
    print separator
    print "Running poszi ..."
    print "This script extract forces from CP2K output file"

    args = main_parser()

    input_file=args.input_file
    output_file=args.output_file

    with open(input_file) as origin_file:
        for line in origin_file:
            f = re.findall(r'\- Atoms:', line)
            if f :
                natm=int(line.split()[-1])
    print "number of atoms : ",natm

    f = open(output_file,'w+')
    ls=None
    with open(input_file) as origin_file:
        for k,line in enumerate(origin_file):
            if len(line.split()) >= 1 : 
                if line.split()[0] =="ATOMIC" and line.split()[1] =="FORCES":
                    ls=k
            if ls :
                if k > ls + 1 and k <= ls + 2 + natm :
                    form = 3*"%6s"+3*"%20s"
                    print >> f , form%(line.split()[0],line.split()[1],line.split()[2],line.split()[3],line.split()[4],line.split()[5], ) 



    print "file "+output_file+" generated"

