#!/usr/bin/python
# ===============================================================
# date       : 30 novembre 2016 
# author     : fmv
# decription : convert POSCAR (vasp) to POSFF (mdff)
# ===============================================================

import datetime
import argparse
from config import Config
from config import Ion

def main_parser():

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i','--input', help='Input file name',required=True,dest="input_file")
    parser.add_argument('-o','--output',help='Output file name',required=True,dest="output_file")

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    separator=60*"="
    now = datetime.datetime.now()

    print(separator)
    print(now.strftime("%Y-%m-%d %H:%M"))
    print("author : filipe.manuel.vasconcelos@gmail.com")
    print(separator)
    print("Running cartoff ...")
    print("This script generate POSFF file from POSCAR file")

    args = main_parser()

    input_file=args.input_file
    output_file=args.output_file

    ions=[]
    conf_pos=Config(ions=ions)

    conf_pos.read_POSCAR(input_file)
    conf_pos.write_POSFF(output_file)
    print("file "+output_file+" generated")


