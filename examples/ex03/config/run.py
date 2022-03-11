#!/usr/bin/python
import sys
sys.path.append('/home/filipe/dev/mdff/bin')
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


    min_energy=float('inf')
    min_it=0
    coeff=0.1
    n=1000
    #2   -1.000000
    #3   -3.000000
    #4   -6.000000
    #5   -9.103852
    #6   -12.712062
    #7   -16.505384
    #8   -19.821489
    #9   -24.113360
    #10   -28.422532
    #11   -32.765970
    #12   -37.967600
    #13   -44.326801
    #natm = 13 
    #target= -44.326801
    natm = 12
    target= -37.967600
    energy = float('inf')

    a = natm ** (1./3.) * 10.0
    print("cluster with",natm,"LJ atoms") 
    print("target energy (M.R. Hoare and P. Pal, Adv. Phys. 20 161 (1971)",target)

    it=0
    while (abs(energy-target)>1e-4) :
        ions=[]
        for ia in range ( natm ) :
            x=(2.0*random.random()-1.0)*a
            y=(2.0*random.random()-1.0)*a
            z=(2.0*random.random()-1.0)*a
            ions.append([x,y,z])
        write_TRAJFF(natm,ions,'TRAJFF')

        cmd = "mdff.x control_opt.F"
        mdff = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        while True:
            line = mdff.stdout.readline().decode('utf-8')
            if line != '':
                if re.search("final energy&pressure =",line.rstrip()):
                    energy = float ( line.rstrip().split()[3] )
                    if energy < min_energy :
                        min_energy = energy
                        min_it = it
                    form="%10s %10i %10s %20.8f %20s %20.8f %10s %10i"
                    print(f"config : {it} energy :{energy} current : {min_energy} config :{min_it}")
            else :
                break
        it+=1
        mdff.stdout.close()
    
    subprocess.call("cp -p ISCFF"+" ISCFF.min", shell=True)
    subprocess.call("cp -p ISTHFF"+" ISTHFF.min", shell=True)

