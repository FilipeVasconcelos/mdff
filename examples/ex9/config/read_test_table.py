#!/usr/bin/python
# ===============================================================
# name        : read_test_table
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
    
    f=open('test.config')

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
        if storing_table and line[0] == "|" and line[1] != "TEST" :
            table.append(line)

        if storing_table : continue

        variables.append(line)

    print "found ",len ( variables ), " variables defined in test.config :"
    print
    for item in variables:
        print item[0]
    print

    print "found ",len ( table), " test defined in test.config :"
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
        if k != len(var) :
            if first :
                print >> fout, "sed -e s/__"+item[0]+"__/"+item[1]+"/g",
            else:
                print >> fout, "    -e s/__"+item[0]+"__/"+item[1]+"/g",
        else:
            print >> fout, "    -e s/__"+item[0]+"__/"+item[1]+"/g tmp1 > tmp"
        first=False
        k+=1

    k=0
    fdo=open('DOALLC','w')
    for item in table:
        k+=1
        fin=open('scr')
        fout=open('scr'+str(k),'w')
        fout.write(fin.read())
        fin.close()
        print "test ",k, item 
    
        print >> fdo, 'Docontrol '+item[0]+' '+item[1]+' '+item[2]+' '+item[3]+' '+item[4]+' '+item[5]


        if item[1] == "YES" : 
            VAL1 = ".true."
        else:
            VAL1=".false."
        if item[2] == "YES" : 
            VAL2 = ".true."
        else:
            VAL2=".false."
        if item[3] == "YES" : 
            VAL3 = ".true."
        else:
            VAL3=".false."

        print >> fout,"sed -e s/__"+"BMH"+"__/"+VAL1+"/g",
        print >> fout,"-e s/__"+"FT_ON"+"__/"+VAL2+"/g",
        print >> fout,"-e s/__"+"COUL"+"__/"+VAL3+"/g tmp",
        fout.close()

    fdo.close()

    return
# ============================================================================
    
var=[]
table=[]
var,table = read_quench()
write_script(var,table)



