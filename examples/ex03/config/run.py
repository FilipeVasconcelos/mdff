#!/usr/bin/python
import sys
sys.path.append('/home/filipe/dev/mdff/bin')
import argparse
import numpy as np
import random
from config import Config
from config import Ion
import subprocess
import re

def write_TRAJFF(natm,ions,filename):

    f=open(filename,'w+')
    f.write(str(natm)+'\n')
    f.write("cluster\n")
    f.write("    500.00000000      0.00000000      0.00000000\n")
    f.write("      0.00000000    500.00000000      0.00000000\n")
    f.write("      0.00000000      0.00000000    500.00000000\n")
    f.write("1\n")
    f.write("A\n")
    f.write(str(natm)+'\n')
    f.write("C\n")
    for ia in range ( natm ) :
        f.write("A "+str(ions[ia][0])+" "+str(ions[ia][1])+" "+str(ions[ia][2])+'\n')
    f.close()
    return

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='Running optimization calculation on LJ cluster')
    parser.add_argument('-n','--natm', help='Number of atoms', required=True,default=12,type=int)
    args = vars(parser.parse_args())
    natm=args["natm"]
    print("cluster with",natm,"LJ atoms") 
    if natm > 13 or natm <3 :
        print("natm should be between 3 and 13 to use target energy")
        target=0.0
    else:
        print("target energy from : M.R. Hoare and P. Pal, Adv. Phys. 20 161 (1971) :",target)
        targets=[None,None,-1.000000 ,-3.000000 ,-6.000000 ,-9.103852 ,\
                           -12.712062,-16.505384,-19.821489,-24.113360,\
                           -28.422532,-32.765970,-37.967600,-44.326801]
        target=targets[natm]

    min_energy=float('inf')
    min_it=0
    nconf=1000
    energy = float('inf')

    a = natm ** (1./3.) * 5.0

    it=0
    while (abs(energy-target)>1e-4) and it < nconf :
        ions=[]
        for ia in range ( natm ) :
            x=(2.0*random.random()-1.0)*a
            y=(2.0*random.random()-1.0)*a
            z=(2.0*random.random()-1.0)*a
            ions.append([x,y,z])
        write_TRAJFF(natm,ions,'TRAJFF')

        cmd = "mpirun -np 2 mdff.x control_opt.F"
        mdff = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        while True:
            line = mdff.stdout.readline().decode('utf-8')
            if line != '':
                if re.search("minimum reached",line.rstrip()) :
                    iters=int(line.rstrip().split()[3])
                if re.search("final energy&pressure =",line.rstrip()):
                    energy = float ( line.rstrip().split()[3] )
                    if energy < min_energy :
                        min_energy = energy
                        min_it = it
                        subprocess.call("cp -p ISCFF"+" ISCFF.min", shell=True)
                        subprocess.call("cp -p ISTHFF"+" ISTHFF.min", shell=True)
                    out=(
                                f'config : {it:4d} '
                                f'energy : {energy:15.8f} ' 
                                f'current (min) : {min_energy:15.8f} '
                                f'iteration : {iters:4d} '
                                f'target : {target:15.8f} '
                        )
                    print(out)
            else :
                break
        it+=1
        mdff.stdout.close()
    print("found target configuration")
    

