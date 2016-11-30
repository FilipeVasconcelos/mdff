#!/usr/bin/env python
# ===============================================================
# date       : 30 nov 2016 
# author     : fmv
# description: This script generate gnuplot plots from input OSZIFF file
# ===============================================================

import matplotlib.pyplot as plt
import numpy as np
import datetime
from optparse import OptionParser


name_quant=['step','Time','Etot','Ekin','Utot','U_vdw','U_coul','Temp','Press','Pvir_vdw','Pvir_coul','Volume','Htot']
units={}
units['Time']='[ps]'
units['Etot']='[eV]'
units['Ekin']='[eV]'
units['Utot']='[eV]'
units['U_vdw']='[eV]'
units['U_coul']='[eV]'
units['Temp']='[K]'
units['Press']='[GPa]'
units['Volume']='[${\AA}^3$]'
units['Htot']='[eV]'

def format_filename(filename):

    f=open(filename)
    first_line = f.readline().split()
    f.close()

    if first_line[0] != "step" and first_line[1] != "=" and first_line[3] != "Time":
        return False
    else:
        return True

def read_OSZIFF(filename):

    alldata=[]
    step=[]
    time=[]
    etot=[]
    ekin=[]
    utot=[]
    uvdw=[]
    ucou=[]
    temp=[]
    pres=[]
    pvir_vdw=[]
    pvir_coul=[]
    volu=[]
    htot=[]

    f=open(filename)
    for i,line in enumerate(f) :
        l = line.split()
        if i%2 == 0:
            step.append(int(l[2]))
            time.append(float(l[5]))
            etot.append(float(l[8]))
            ekin.append(float(l[11]))
            utot.append(float(l[14]))
            uvdw.append(float(l[17]))
            ucou.append(float(l[20]))
        if i%2 == 1:
            temp.append(float(l[8]))
            pres.append(float(l[11]))
            pvir_vdw.append(float(l[14]))
            pvir_coul.append(float(l[17]))
            volu.append(float(l[20]))
            htot.append(float(l[23]))

    f.close()

    alldata.append(step)
    alldata.append(time)
    alldata.append(etot)
    alldata.append(ekin)
    alldata.append(utot)
    alldata.append(uvdw)
    alldata.append(ucou)
    alldata.append(temp)
    alldata.append(pres)
    alldata.append(pvir_vdw)
    alldata.append(pvir_coul)
    alldata.append(volu)
    alldata.append(htot)

    return alldata
    
def averaging(name_quant,alldata,last_points):

    average=[]
    average.append(None)#step
    average.append(None)#time
    for i,l in enumerate(alldata):
        if name_quant[i] != "step" and name_quant[i] != "Time":
            average.append(np.mean(l[:last_points]))
            print '{0:<11} {1:<10} {2:15.8e} {3:^10} {4:15.8e} '.format("<"+name_quant[i]+">","=",np.mean(l[:last_points]),"std.",np.std(l[:last_points]))

    return average

def plot_quant2(alldata,name_quant,average,q,title):

    t = alldata[1] 
    s = alldata[q[0]] 
    x1= [average[q[0]]]*len(t)
    plt.plot(t, s,'b')
    plt.plot(t, x1, '--' , color='b') 

    s = alldata[q[1]] 
    x1= [average[q[1]]]*len(t)
    plt.plot(t, s , 'g')
    plt.plot(t, x1, '--', color='g') 

    plt.xlabel('Time'+' '+units['Time'])
    plt.ylabel(name_quant[q[0]]+' '+units[name_quant[q[0]]])
    plt.title(title)
    plt.show()

def plot_quant(alldata,name_quant,average,q,title):

    t = alldata[1] 
    s = alldata[q] 
    x1= [average[q]]*len(t)
    plt.plot(t, s,'b')
    plt.plot(t, x1, '--' , color='b') 

    plt.xlabel('Time'+' '+units['Time'])
    plt.ylabel(name_quant[q]+' '+units[name_quant[q]])
    plt.title(title)
    plt.show()


def main_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",version="%prog 0.1")
    parser.add_option("-n", "--no_plot",dest="plot_flag",action="store_false",default=True,
            help="show plot of average and instantaneous thermodynamic quantities along OSZIFF")
    parser.add_option("-l", "--last",dest="last_points",default=None,
            help="averaging on last <l> points")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("wrong number of arguments")

    return options,args

if __name__ == '__main__':

    separator=60*"="
    now = datetime.datetime.now()

    print separator
    print now.strftime("%Y-%m-%d %H:%M")
    print "author : filipe.manuel.vasconcelos@gmail.com"
    print separator
    print "Running poszi ..."
    print "This script generate gnuplot plots from input OSZIFF file"

    options,args = main_parser()
    filename=args[0]
    option_dict = vars(options)


    if format_filename(filename):
        alldata = read_OSZIFF(filename)
    else:
        raise ValueError('Format of input file doesnt match OSZIFF')

    last_points = option_dict['last_points'] 
    plot_flag   = option_dict['plot_flag'] 

    if last_points == None :
        last_points = len(alldata[0])
    else:
        last_points = int( last_points ) 

    print len(alldata[0])," points in input file"
    print "averaging on last",last_points,"points"

    average=averaging(name_quant,alldata,last_points)
   
    if plot_flag :
        plot_quant2(alldata,name_quant,average,[2,12],title='Total Energy MD')
        plot_quant(alldata,name_quant,average,4,title='Potential Energy MD')
        plot_quant(alldata,name_quant,average,7,title='Temperature MD')
        plot_quant(alldata,name_quant,average,11,title='Volume MD')
