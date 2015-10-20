#!/usr/bin/python
# ===============================================================
# name        : read_quench_table
# author      : fmv
# date        : somewhere in march-april 2014
# description : this python script reads the quench.config global
#               parameters and those related to the quench. Then 
#               it generates "sed" scripts according to the input
#               parameters together with the DOALLC DOALLMD exec-
#               ution sequence.
# ===============================================================


# ------------------------------------------------
# main function to read the quench.config file
# it returns variables , newtable which are list of 
# the main lines of quench.config.
# variables : global parameters 
#             for item in var:
#               item[0] = name
#               item[1] = value
# newtable  : quench table
# ------------------------------------------------
def read_quench():
    
    f=open('quench.config')

    lines=f.readlines()
    f.close()

    storing_table = False
    table=[]
    variables=[]
    for item in lines:
        line=[]
        line=item.split()
        if len(line) < 1 : continue
        if line[0] == "TABLE:" :
            storing_table = True
        if storing_table and line[0] == "|" and line[1] != "TEMP" :
            table.append(line)

        if storing_table : continue

        variables.append(line)

    print "found ",len ( variables ), " variables defined in quench.config :"
    print
    for item in variables:
        print item[0]
    print

#+--------+--------+--------------+-------------+----------+--------------+------+------------+
#|  TEMP  |   ENS. |  npas [ps]   | fprint [ps] |  traj    |  period [ps] |  POT |    flag    |
#+--------+--------+--------------+-------------+----------+--------------+------+------------+

    print "found ",len ( table), " runs defined in quench.config :"
    newtable=[]
    for item in table:
        newline=[]
        newline=filter(lambda x: x != "|", item)
        newtable.append(newline)

    return variables,newtable

# ------------------------------------------------
# ------------------------------------------------
def write_script(var,table):

    # first write the sed script for the global parameters 
    # no matters the equivalent parameter exists in the template 
    fout=open('scr','w')
    first=True
    k=1
    for item in var:
        if item[0] == "DT" : dt=float(item[1])
        if k != len(var) :
            if first :
                print >> fout, "sed -e s/__"+item[0]+"__/"+item[1]+"/g",
            else:
                print >> fout, "    -e s/__"+item[0]+"__/"+item[1]+"/g",
        else:
            print >> fout, "    -e s/__"+item[0]+"__/"+item[1]+"/g tmp1 > tmp"
        first=False
        k+=1
    fout.close()


    fdo=open('DOALLC','w')
    fdo2=open('DOALLMD','w')
    first=True
    k=0
    for item in table:

        k+=1
        fin=open('scr')
        fout=open('scr'+str(k),'w')
        fout.write(fin.read())
        fin.close()
        print "run ",k, item 
        if item[0] != "QUENCH" :

            NAME="TEMP"
            valtemp=item[0] 
            VAL=str(float(valtemp[:-1]))
            print >> fout, "sed    -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="NPAS"
            VAL=str(int(float(item[2])/dt))
            print >> fout,  "    -e s/__"+NAME+"__/"+VAL+"/g",
    
            NAME="LTRAJ"
            VAL=item[4]
            print >> fout,  "   -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="INTEGRATOR"
            valtemp=item[1]
            if valtemp == "NVE"   : VAL="\\'nve-vv\\'"
            if valtemp == "NVT"   : VAL="\\'nvt-nhcn\\'"
            if valtemp == "NPT_I" : VAL="\\'npt-nhcnp\\'"
            print >> fout,  "   -e s/__"+NAME+"__/"+VAL+"/g",

            if item[1] == "NVE" and item[7] == "rescaling" :
                NAME="NEQUIL" 
                VAL=str(int(float(item[2])/dt*0.9))
                print >> fout,  "   -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="ANNEALING"
            VAL="1.0000"
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",
           
            NAME="NPRINT"
            VAL=str(int(float(item[3])/dt))
            print >> fout,  "   -e s/__"+NAME+"__/"+VAL+"/g",
    
            NAME="FPRINT"
            VAL=str(int(float(item[3])/dt))
            print >> fout,  "   -e s/__"+NAME+"__/"+VAL+"/g",
    
            NAME="ITRAJ_PERIOD"
            if item[5] != "-":
                VAL=str(int(float(item[5])/dt))
            else:
                VAL=str(int(float(item[2])/dt))
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

            if first :
                NAME="RESTART"
                VAL=".false."
                print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

                NAME="DATA"
                VAL="\\'rnn\\'"
                print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g tmp " 
            else:
                NAME="RESTART"
                VAL=".true."
                print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

                NAME="DATA"
                VAL="\\'rvf\\'"
                print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g tmp"
    
    
            print >> fdo, 'Docontrol "'+item[1]+' '+item[6]+' '+item[0]+'" control'+'_'+item[1]+'_'+item[6]+'_'+item[0]+'_'+item[7]+'.F'\
                          +' '+item[1]+' '+item[6]+' '+item[7]+' '+str(k)
            print >> fdo2, 'DoMD "'+item[1]+' '+item[6]+' '+item[0]+'" control'+'_'+item[1]+'_'+item[6]+'_'+item[0]+'_'+item[7]+'.F'\
                          +' '+item[1]+' '+item[6]+' '+item[7]+' '+str(k)+' '+item[1]+'_'+item[0]+'_'+item[6]+'_'+item[7]
        # QUENCH
        else:

            NAME="TEMP"
            VAL=""
            print >> fout, "sed -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="NPAS"
            VAL=str(int(float(item[2])/dt))
            print >> fout, "    -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="LTRAJ"
            VAL=item[4]
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="ITRAJ_PERIOD"
            if item[5] != "-":
                VAL=str(int(float(item[5])/dt))
            else:
                VAL=str(int(float(item[2])/dt))
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",


            NAME="INTEGRATOR"
            VAL="\\'nve-vv\\'"
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="NPRINT"
            VAL=str(int(float(item[3])/dt))
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",
    
            NAME="FPRINT"
            VAL=str(int(float(item[3])/dt))
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",
               
            NAME="ANNEALING"
            VAL=item[7]
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",
    
            NAME="RESTART"
            VAL=".true."
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g",

            NAME="DATA"
            VAL="'rvf'"
            print >> fout, "   -e s/__"+NAME+"__/"+VAL+"/g tmp"


            print >> fdo, 'Docontrol "'+item[1]+' '+item[6]+' '+item[0]+'" control'+'_'+item[1]+'_'+item[6]+'_'+item[0]+'.F'+' '+item[1]+' '+item[6]+' '+item[0]+' '+str(k) 
            print >> fdo2,'DoMD "'     +item[1]+' '+item[6]+' '+item[0]+'" control'+'_'+item[1]+'_'+item[6]+'_'+item[0]+'.F'+' '+item[1]+' '+item[6]+' '+item[0]+' '+str(k)+' '+item[1]+'_'+item[0]+'_'+item[6]+'_'+item[6]
        first=False

    fout.close()
    fdo.close()
    fdo2.close()

    return
# ============================================================================
    
var=[]
table=[]
var,table = read_quench()
write_script(var,table)



