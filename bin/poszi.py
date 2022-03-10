#!/usr/bin/env python3
# ===============================================================
# date       : mars 2022
# author     : fmv
# description: This script generate matplotlib plots from input OSZIFF file
# ===============================================================
import matplotlib.pyplot as plt
import numpy as np
import datetime
import argparse
import sys

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
# ------------------------------------------------------------------------------
def format_filename(filename):

    f=open(filename)
    first_line = f.readline().split()
    f.close()
    if first_line[0] != 110*"-":
        return False
    else:
        return True
# ------------------------------------------------------------------------------
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

    etot_l=open("etot_l","w+") 
    temp_l=open("temp_l","w+") 
    f=open(filename)
    k=0
    for i,line in enumerate(f) :
        l = line.split()
        if l[0] == 110*'-' : continue
        if l[0] == 'step' : continue
        if k%2 == 0 :
            step.append(int(l[0]))
            time.append(float(l[1]))
            etot.append(float(l[2]))
            ekin.append(float(l[3]))
            utot.append(float(l[4]))
            uvdw.append(float(l[5]))
            ucou.append(float(l[6]))
            etot_l.write(line)
        if k%2 == 1:
            temp.append(float(l[1]))
            pres.append(float(l[2]))
            pvir_vdw.append(float(l[3]))
            pvir_coul.append(float(l[4]))
            volu.append(float(l[5]))
            htot.append(float(l[6]))
            temp_l.write(line)
        k+=1
    f.close()
    etot_l.close()
    temp_l.close()
    if len(step) > len(temp) :
        print("last line is missing")
        step.pop()
        time.pop()
        etot.pop()
        ekin.pop()
        utot.pop()
        uvdw.pop()
        ucou.pop()
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
# ------------------------------------------------------------------------------
def averaging(name_quant,alldata,last_points):
    c='>'
    average=[]
    average.append(None)#step
    average.append(None)#time
    print(50*'-')
    for i,l in enumerate(alldata):
        ln=len(name_quant[i])
        if name_quant[i] != "step" and name_quant[i] != "Time":
            average.append(np.mean(l[-last_points:]))
            print(f"<{name_quant[i]:<s}{c:{10-ln}s} = {np.mean(l[-last_points:]):15.8e} std. {np.std(l[-last_points:]):<15.8e}")
    print(50*'-')
    return average

def plot_data(alldata,name_quant,average,last_points):

    data_to_plot=[{"title":"Total Energy",
                   "index": [2,12],
                   "color":["blue","green"],
                   "grid": [0,0]},
                  {"title":"Potential Energy (total)",
                   "index": 4,
                   "color":"blue",
                   "grid": [1,0]},
                  {"title":"Potential Energy (Uvdw)",
                   "index": 5,
                   "color":"blue",
                   "grid": [1,1]},
                  {"title":"Potential Energy (Ucoul)",
                   "index": 6,
                   "color":"blue",
                   "grid": [1,2]},
                  {"title":"Temperature",
                   "index": 7,
                   "color": "blue",
                   "grid": [0,1]},
                  {"title":"Volume",
                   "index": 11,
                   "color": "blue",
                   "grid": [0,2]}
                 ]

    fig, axs = plt.subplots(2, 3,figsize=(12,6))
    t=alldata[1][-last_points:]
    for data in data_to_plot:
        spi,spj=data["grid"][0],data["grid"][1]
        title=data["title"]
        #print(spi,spj,title)
        if (isinstance(data["index"],list)) :
            for k,index in enumerate(data["index"]):
                color=data["color"][k]
                s=alldata[index][-last_points:]
                mean=[average[index]]*len(t)
                axs[spi][spj].plot(t,s,color)
                axs[spi][spj].plot(t,mean,'--',color=color)
                axs[spi][spj].set_xlabel('Time'+' '+units['Time'])
                axs[spi][spj].set_ylabel(name_quant[index]+' '+units[name_quant[index]])
                axs[spi][spj].set_title(title)
        else:
            index=data["index"]
            color=data["color"]
            s=alldata[index][-last_points:]
            mean=[average[index]]*len(t)
            axs[spi][spj].plot(t,s,color)
            axs[spi][spj].plot(t,mean,'--',color=color)
            axs[spi][spj].set_xlabel('Time'+' '+units['Time'])
            axs[spi][spj].set_ylabel(name_quant[index]+' '+units[name_quant[index]])
            axs[spi][spj].set_title(title)

    fig.tight_layout()
    plt.show()

def main_parser():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i','--input', help='Input file name',required=True,dest="input_file")
    parser.add_argument("-n", "--no_plot",dest="plot_flag",action="store_false",default=True,
            help="show plot of average and instantaneous thermodynamic quantities along OSZIFF")
    parser.add_argument("-l", "--last",dest="last_points",default=None,
            help="averaging on last <l> points")

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    separator=60*"="
    now = datetime.datetime.now()

    print(separator)
    print(now.strftime("%Y-%m-%d %H:%M"))
    print("author : filipe.manuel.vasconcelos@gmail.com")
    print(separator)
    print("Running poszi ...")
    print("This script generate matplotlib plots from input OSZIFF file")

    args = main_parser()
    input_file=args.input_file

    if format_filename(input_file):
        alldata = read_OSZIFF(input_file)
    else:
        raise ValueError('Format of input file doesnt match OSZIFF')

    if args.last_points == None :
        last_points = len(alldata[0])
    else:
        last_points = int( args.last_points ) 

    print(f"{len(alldata[0])} points in input file")
    if len(alldata[0]) != 0 : 
        average=averaging(name_quant,alldata,last_points)
        print("averaging",end=' ')
        if args.plot_flag :
            print("and plot",end=' ')
            plot_data(alldata,name_quant,average,last_points)
        print(f"on last {last_points} points")

