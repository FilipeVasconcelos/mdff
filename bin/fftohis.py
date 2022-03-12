#!/usr/bin/python3
# ===============================================================
# date       : 11 mars 2022 
# author     : fmv
# decription : convert TRAJFF (mdff) to HISTORY (DLPOLY)
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
    print("Running fftohis ...")
    print("This script generate HISTORY file (DLPOLY) from TRAJFF file (MDFF)")

    args = main_parser()

    input_file=args.input_file
    output_file=args.output_file

    ions=[]
    conf_pos=Config(ions=ions)

    conf_pos.read_TRAJFF_write(input_file,output_file)
    print("file "+output_file+" generated")


