#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
# 
import Tkinter as tk
import ttk 
import tkFileDialog
import tkMessageBox
import tkFont
import sys
import os
import signal
import platform
import subprocess
from multilistbox import MultiListbox
import matplotlib
import random
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from distutils.spawn import find_executable

# local libs
from create_control_file import get_control_file, namelist_all, all_values_allowed, integrator_allowed,default_values
from poszi import read_OSZIFF,averaging,format_filename,plot_quant2,plot_quant,units
from config import Config
from config import Ion

class MApp(tk.Frame):

    #----------------------------------------------------------------------
    def __init__(self, parent):
        """Constructor"""
        tk.Frame.__init__(self,parent)
        self.parent = parent
        self.parent.title('MDff GUI'+"~")
        w, h = self.parent.winfo_screenwidth(), root.winfo_screenheight()
        self.parent.geometry("%dx%d+0+0" % (w, h))
        self.pack()
        self.global_path='/home/filipe/dev/github/mdff'
        self.doc_path=self.global_path+'/doc'
        self.exe_path=tk.StringVar(None)
        self.exe_path.set(find_executable('mdff.x'))
        self.control_path=tk.StringVar(None)
        self.control_path.set('./control.F')
        self.posff_path=tk.StringVar(None)
        self.posff_path.set('./POSFF')
        self.stdout_path=tk.StringVar(None)
        self.stdout_path.set('./log')

        self.running=True
    
        self._create_panel_files()
        self._create_panel_run()
        self._create_menu()
        # shortscut
        self.parent.bind_all("<Control-w>", self._quit)
        self.parent.bind_all("<Control-s>", self._save)
        self.parent.bind_all("<Control-x>", self.start_exe)
        self.parent.bind_all("<Control-k>", self._kill_exe)


        self.run_frame = None
        self.plot_frame = None
        self.config_frame = None
        self.control_md_frame = None

    #----------------------------------------------------------------------
    def _create_menu(self):

        self.menubar = tk.Menu(self)

        menuFile = tk.Menu(self,tearoff=0)

        menuCreate = tk.Menu(menuFile,tearoff=0)

        menuCreate.add_command(label="full",command=lambda :self._create_control_file("full"))
        menuCreate.add_command(label="empty",command=lambda :self._create_control_file("empty"))
        menuCreate.add_command(label="md",command=lambda :self._create_control_file("md"))
        menuCreate.add_command(label="opt",command=lambda :self._create_control_file("opt"))
        menuCreate.add_command(label="stochio",command=lambda :self._create_control_file("stochio"))

        menuFile.add_cascade(label="New control file",menu=menuCreate,underline=0)
        menuFile.add_command(label="New POSFF file",command=self._create_posff_file)
        menuFile.add_command(label="Save as",command=self._save_as)
        menuFile.add_separator()
        menuFile.add_command(label="Quit",command=self.parent.destroy, accelerator="Ctrl+w")

        self.menubar.add_cascade(label="File", menu=menuFile)

        menuExe = tk.Menu(self,tearoff=0)
        menuExe.add_command(label="Select control file", command=self._select_control_file)
        menuExe.add_command(label="Select POSFF file", command=self._select_posff_file)
        menuExe.add_command(label="Select executable",command=self._select_exe)
        menuExe.add_command(label="Execute",command=self._start_exe, accelerator="Ctrl+x")

        self.menubar.add_cascade(label="Run", menu=menuExe)

        menuAnalysis = tk.Menu(self,tearoff=0) 
        menuAnalysis.add_command(label="Thermodynamic plots", command=self._run_stat)
        menuAnalysis.add_command(label="Radial Distribution Function", command=self._select_posff_file)
        self.menubar.add_cascade(label="Analysis", menu=menuAnalysis)

        menuHelp = tk.Menu(self,tearoff=0)
        menuHelp.add_command(label="Documentation",command=self._doc)
        menuHelp.add_command(label="About",command=self._apropos)
        self.menubar.add_cascade(label="Help", menu=menuHelp)


        try:
            self.parent.config(menu=self.menubar)
        except AttributeError:
            self.parent.call(parent, "config", "-menu", self.menubar)

    #----------------------------------------------------------------------
    def _create_panel_run(self):
        self.run_frame = tk.Frame(self.parent, name='run_frame')
        self.run_frame.pack(side=tk.TOP,fill=tk.X)

        lbl=tk.LabelFrame(self.run_frame,text="Input general information ", container=False,width=1000)
#        lbl.grid(row=0, column=0,sticky=tk.W+tk.E,padx=5,pady=10)
        lbl.pack(side=tk.TOP,fill=tk.X,padx=10,pady=10)
        # control file 
        lc = tk.Label(lbl, text="current control file : " ,justify = tk.LEFT)
        lc.grid(row=0,column=0,sticky=tk.W,padx=5)
        mc=tk.Message(lbl, textvariable=self.control_path , justify=tk.LEFT,width=1000)
        mc.grid(row=0,column=1,sticky=tk.W,pady=1,padx=5)
        # POSFF file 
        lp = tk.Label(lbl, text="current config file : ",justify = tk.LEFT)
        lp.grid(row=1,column=0,sticky=tk.W,padx=5)
        mp=tk.Message(lbl, textvariable=self.posff_path , justify=tk.LEFT,width=1000)
        mp.grid(row=1,column=1,sticky=tk.W,pady=1,padx=5)
        # stdout file 
        ls = tk.Label(lbl, text="standard output (std): ",justify = tk.LEFT)
        ls.grid(row=2,column=0,sticky=tk.W,padx=5)
        ms=tk.Message(lbl, textvariable=self.stdout_path , justify=tk.LEFT,width=1000)
        ms.grid(row=2,column=1,sticky=tk.W,pady=1,padx=5)
        # exe file 
        lx = tk.Label(lbl, text="executable path : ",justify = tk.LEFT)
        lx.grid(row=3,column=0,sticky=tk.W,padx=5)
        mx=tk.Message(lbl, textvariable=self.exe_path , justify=tk.LEFT,width=1000)
        mx.grid(row=3,column=1,sticky=tk.W,pady=1,padx=5)



    #----------------------------------------------------------------------
    def _create_panel_files(self):
        self.files_frame = tk.Frame(self.parent, name='files_frame')
        self.files_frame.pack(side=tk.LEFT,fill=tk.Y)

        # create the notebook
        self.nb = ttk.Notebook(self.files_frame, name='notebook')
        # extend bindings to top level window allowing
        #   CTRL+TAB - cycles thru tabs
        #   SHIFT+CTRL+TAB - previous tab
        #   ALT+K - select tab using mnemonic (K = underlined letter)
        self.nb.enable_traversal()
        self.nb.pack(side=tk.LEFT,fill=tk.Y)
        self._create_control_tab()
        self._create_posff_tab()
        self._create_stdout_tab()

    #----------------------------------------------------------------------
    def _create_stdout_tab(self):
        self.frame_stdout = tk.Frame(self.nb, name='log')
        self.frame_stdout.pack(side=tk.LEFT,fill=tk.Y)
        self.nb.add(self.frame_stdout, text='log', underline=0)
        
        self.txt_stdout = tk.Text(self.frame_stdout,width=100,height=44)#,background="gray50") 
        vscroll = ttk.Scrollbar(self.frame_stdout, orient=tk.VERTICAL, command=self.txt_stdout.yview)
        self.txt_stdout['yscroll'] = vscroll.set
        self.txt_stdout.grid(row=0,column=0,sticky=tk.W)
        vscroll.grid(row=0, column=1,sticky=tk.N+tk.S,padx=5)
        

    #----------------------------------------------------------------------
    def _create_control_tab(self):

        self.frame_control = tk.Frame(self.nb, name='control')
        self.frame_control.pack(side=tk.LEFT,fill=tk.Y)
        self.nb.add(self.frame_control, text='control.F', underline=0)

        self.txt_control = tk.Text(self.frame_control,width=100,height=44,undo=True)#,background="gray50") 
        vscroll = ttk.Scrollbar(self.frame_control, orient=tk.VERTICAL, command=self.txt_control.yview)
        self.txt_control['yscroll'] = vscroll.set
        self.txt_control.grid(row=0, column=0,sticky=tk.N+tk.S)
        vscroll.grid(row=0, column=1,sticky=tk.N+tk.S,padx=5)


        self.txt_control.tag_configure("tag", foreground="navy")
        self.txt_control.tag_configure("int", foreground="red")
        self.txt_control.tag_configure("tag2", foreground='green4')
        self.txt_control.tag_configure("comment", foreground="black")
        self.txt_control.bind("<Any-KeyRelease>", self.highlight_control)
        self.txt_control.bind("<Any-ButtonRelease>", self.highlight_control)

        file = open(self.control_path.get(), 'a+')
        if file != None:
            content = file.read()
            self.txt_control.delete('1.0','end')
            self.txt_control.insert('1.0',content)
            file.close()
        self.txt_control.update()

    #----------------------------------------------------------------------
    def highlight_control(self, event=None):

        '''
            put some colors on control panel
        '''

        #tags highlight 
        counts = []
        self.txt_control.tag_remove("tag", "1.0", "end")
        keys = ["controltag","mdtag","fieldtag","stochiotag","end"]

        for key in keys:
            counts.append(tk.IntVar())
            self.txt_control.mark_set("matchStart", "1.0")
            self.txt_control.mark_set("matchEnd", "1.0")
            while True:
                index = self.txt_control.search(key, "matchEnd","end", count=counts[-1])
                if index == "" : break # no match was found

                self.txt_control.mark_set("matchStart", index)
                self.txt_control.mark_set("matchEnd", "%s+%sc" % (index, counts[-1].get()))
                self.txt_control.tag_add("tag", "matchStart", "matchEnd")

        #values (numbers, boolean)
        counts = []
        self.txt_control.tag_remove("int", "1.0", "end")
        keys = ["1","2",'3','4','5','6','7','8','9','0','d0','_dp','.true.','.false.'] + all_values_allowed

        for key in keys:
            counts.append(tk.IntVar())
            self.txt_control.mark_set("matchStart", "1.0")
            self.txt_control.mark_set("matchEnd", "1.0")
            while True:
                index = self.txt_control.search(key, "matchEnd","end", count=counts[-1])
                if index == "" : break # no match was found

                self.txt_control.mark_set("matchStart", index)
                self.txt_control.mark_set("matchEnd", "%s+%sc" % (index, counts[-1].get()))
                self.txt_control.tag_add("int", "matchStart", "matchEnd")

        # parameters name 
        counts = []
        self.txt_control.tag_remove("tag2", "1.0", "end")
        keys = namelist_all

        for key in keys:
            counts.append(tk.IntVar())
            self.txt_control.mark_set("matchStart", "1.0")
            self.txt_control.mark_set("matchEnd", "1.0")
            while True:
                index = self.txt_control.search(key, "matchEnd","end", count=counts[-1])
                if index == "" : break # no match was found

                self.txt_control.mark_set("matchStart", index)
                self.txt_control.mark_set("matchEnd", "%s+%sc" % (index, counts[-1].get()))
                self.txt_control.tag_add("tag2", "matchStart", "matchEnd")

        #comments
        self.txt_control.tag_remove("comment", "1.0", tk.END)

        counts=tk.IntVar()
        self.txt_control.mark_set("matchStart", "1.0")
        self.txt_control.mark_set("matchEnd", "1.0")
        while True:
            index = self.txt_control.search("!", "matchEnd","end", count=counts)
            if index == "" :  break
            self.txt_control.mark_set("matchStart", index)
            self.txt_control.mark_set("matchEnd", "%s %s" % (index, "lineend"))
            self.txt_control.tag_add("comment", "matchStart" , "matchEnd")

    #----------------------------------------------------------------------
    def _create_posff_tab(self):

        self.frame_posff = tk.Frame(self.nb, name='posff')
        self.frame_posff.pack(side=tk.LEFT,fill=tk.Y)
        self.nb.add(self.frame_posff, text='POSFF', underline=0)

        self.txt_posff = tk.Text(self.frame_posff, width=100,height=44,wrap=tk.NONE)#,background="gray50") 
        yscroll = ttk.Scrollbar(self.frame_posff, orient=tk.VERTICAL, command=self.txt_posff.yview)
        xscroll = ttk.Scrollbar(self.frame_posff, orient=tk.HORIZONTAL, command=self.txt_posff.xview)
        self.txt_posff['yscroll'] = yscroll.set
        self.txt_posff['xscroll'] = xscroll.set
        self.txt_posff.grid(row=0,column=0,sticky=tk.W)
        xscroll.grid(row=1,column=0,sticky=tk.W+tk.E,pady=5)
        yscroll.grid(row=0,column=2,sticky=tk.N+tk.S,padx=5)

        file = open(self.posff_path.get(), 'rb')
        if file != None:
            content = file.read()
            self.txt_posff.delete('1.0','end')
            self.txt_posff.insert('1.0',content)
            file.close()


    #----------------------------------------------------------------------
    def _create_control_file(self,calc):
        filename='control_'+calc+'.F'
        if calc in ["full",'empty','stochio'] :
            get_control_file(calc,default_values,None,default_values,None,default_values,None)
            file = open(filename, 'rb')
            if file != None:
                content = file.read()
                self.txt_control.delete('1.0','end')
                self.txt_control.insert('1.0',content)
                file.close()
        elif calc == "md":
            self.md_values = default_values.copy()
            self.field_values = default_values.copy()
            self.control_values = default_values.copy()
            # remove possible frame in the same region
            if self.plot_frame is not None :
                self.plot_frame.destroy()
            if self.config_frame is not None :    
                self.config_frame .destroy()

            self.control_md_frame = tk.Frame(self.parent, name='control_md_frame')
            self.control_md_frame.pack(side=tk.LEFT,fill=tk.BOTH)
#            mainlabel = tk.Label(self.control_md_frame, text="Construct control.F for calc='md' \n(mind only few parameters can be access from this GUI)")
#            mainlabel.grid(row=0,column=0,sticky=tk.W+tk.E,pady=5,padx=10)
       

            # integrator
            lintegrator = tk.Label(self.control_md_frame, text="integrator")
            lintegrator.grid(row=1,column=0,sticky=tk.W+tk.E,pady=5,padx=5)
            self.cv_integrator={'NVE':tk.IntVar(),'NVT':tk.IntVar(),'NPT':tk.IntVar()}
            self.cv_integrator['NVE'].set(0)
            self.cv_integrator['NVT'].set(0)
            self.cv_integrator['NPT'].set(0)

            self.cb_integrator_nve = tk.Checkbutton(self.control_md_frame, text='NVE', variable=self.cv_integrator['NVE'],  onvalue = 1, offvalue = 0 , command=self._cb_integrator_nve)
            self.cb_integrator_nve.grid(row=1,column=1,sticky=tk.W+tk.E,pady=5,padx=5)

            self.cb_integrator_nvt = tk.Checkbutton(self.control_md_frame, text='NVT', variable=self.cv_integrator['NVT'],  onvalue = 1, offvalue = 0 , command=self._cb_integrator_nvt)
            self.cb_integrator_nvt.grid(row=1,column=2,sticky=tk.W+tk.E,pady=5,padx=5)

            self.cb_integrator_npt = tk.Checkbutton(self.control_md_frame, text='NPT', variable=self.cv_integrator['NPT'],  onvalue = 1, offvalue = 0 , command=self._cb_integrator_npt)
            self.cb_integrator_npt.grid(row=1,column=3,sticky=tk.W+tk.E,pady=5,padx=5)


            #potential
            lpot = tk.Label(self.control_md_frame, text="potential")
            lpot.grid(row=2,column=0,sticky=tk.W+tk.E,pady=5,padx=5)
            self.cv_pot={'LJ':tk.IntVar(),'Morse':tk.IntVar(),'RI':tk.IntVar(),'BMHFTD':tk.IntVar(),'PIM':tk.IntVar()}
            self.cv_pot['LJ'].set(0)
            self.cv_pot['Morse'].set(0)
            self.cv_pot['RI'].set(0)
            self.cv_pot['BMHFTD'].set(0)
            self.cv_pot['PIM'].set(0)
    
            self.cb_pot_lj = tk.Checkbutton(self.control_md_frame, text='Lennard-Jones', variable=self.cv_pot['LJ'],  onvalue = 1, offvalue = 0 , command=self._cb_pot_lj)
            self.cb_pot_lj.grid(row=2,column=1,sticky=tk.W+tk.E,pady=5,padx=5)

            self.cb_pot_morse = tk.Checkbutton(self.control_md_frame, text='Morse', variable=self.cv_pot['Morse'],  onvalue = 1, offvalue = 0 , command=self._cb_pot_morse)
            self.cb_pot_morse.grid(row=2,column=2,sticky=tk.W+tk.E,pady=5,padx=5)

            self.cb_pot_ri = tk.Checkbutton(self.control_md_frame, text='Rigid Ion', variable=self.cv_pot['RI'],  onvalue = 1, offvalue = 0 , command=self._cb_pot_ri)
            self.cb_pot_ri.grid(row=2,column=3,sticky=tk.W+tk.E,pady=5,padx=5)


            self.cb_pot_bmhftd = tk.Checkbutton(self.control_md_frame, text='BMHFTD', variable=self.cv_pot['BMHFTD'],  onvalue = 1, offvalue = 0 , command=self._cb_pot_bmhftd)
            self.cb_pot_bmhftd.grid(row=2,column=4,sticky=tk.W+tk.E,pady=5,padx=5)

            self.cb_pot_pim = tk.Checkbutton(self.control_md_frame, text='PIM', variable=self.cv_pot['PIM'],  onvalue = 1, offvalue = 0 , command=self._cb_pot_pim)
            self.cb_pot_pim.grid(row=2,column=5,sticky=tk.W+tk.E,pady=5,padx=5)

            gene=tk.Button(self.control_md_frame, text="Generate control file", command = lambda:self._gene_md(calc))
            gene.grid(row=9,column=0,columnspan=5,sticky=tk.W+tk.E,pady=5)

        
        self.nb.select(0) # show control tab
    #----------------------------------------------------------------------
    def _gene_md(self,calc):
        filename='control_'+calc+'.F'
        kintegration=None
        kcontrol=None
        kpotential=None

        if self.cv_integrator['NVE'].get() == 1 :
            kintegration='NVE'
        if self.cv_integrator['NVT'].get() == 1 :
            kintegration='NVT'
        if self.cv_integrator['NPT'].get() == 1 :
            kintegration='NPT'
        if self.cv_pot['LJ'].get() == 1 :
            kpotential = 'LJ'
            kcontrol='LJ'
        if self.cv_pot['Morse'].get() == 1 :
            kpotential = 'Morse'
            kcontrol='Morse'
        if self.cv_pot['PIM'].get() == 1 :
            kpotential = 'PIM'
            kcontrol='PIM'
        if self.cv_pot['RI'].get() == 1 :
            kpotential = 'RI'
            kcontrol='RI'
        if self.cv_pot['BMHFTD'].get() == 1 :
            kpotential = 'BMHFTD'
            kcontrol='BMHFTD'
        get_control_file(calc,self.control_values,kcontrol,self.md_values,kintegration,self.field_values,kpotential)
        file = open(filename, 'rb')
        if file != None:
            content = file.read()
            self.txt_control.delete('1.0','end')
            self.txt_control.insert('1.0',content)
            file.close()
        self.nb.select(0) # show control tab

    #----------------------------------------------------------------------
    def _cb_pot_lj(self,event=None):
        if self.cv_pot['LJ'].get() == 0:
            self.cb_pot_morse.configure(state="normal")
            self.cb_pot_ri.configure(state="normal")
            self.cb_pot_pim.configure(state="normal")
            self.cb_pot_bmhftd.configure(state="normal")
        else:
            self.cb_pot_morse.configure(state=tk.DISABLED)
            self.cb_pot_ri.configure(state=tk.DISABLED)
            self.cb_pot_pim.configure(state=tk.DISABLED)
            self.cb_pot_bmhftd.configure(state=tk.DISABLED)
            self.control_values['lnmlj'] = ".true."
            self.control_values['lvnlist'] = ".true."
            self.control_values['lreduced'] = ".true."


    #----------------------------------------------------------------------
    def _cb_pot_ri(self,event=None):
        if self.cv_pot['RI'].get() == 0:
            self.cb_pot_morse.configure(state="normal")
            self.cb_pot_lj.configure(state="normal")
            self.cb_pot_pim.configure(state="normal")
            self.cb_pot_bmhftd.configure(state="normal")
        else:
            self.cb_pot_morse.configure(state=tk.DISABLED)
            self.cb_pot_lj.configure(state=tk.DISABLED)
            self.cb_pot_pim.configure(state=tk.DISABLED)
            self.cb_pot_bmhftd.configure(state=tk.DISABLED)
            self.field_values['lautoES'] = ".true."
            self.control_values['lcoulomb'] = ".true."
            self.control_values['lbmhftd'] = ".true."
            self.control_values['lbmhft'] = ".true."
            self.control_values['lreduced'] = ".true."


    #----------------------------------------------------------------------
    def _cb_pot_morse(self,event=None):
        if self.cv_pot['Morse'].get() == 0:
            self.cb_pot_lj.configure(state="normal")
            self.cb_pot_ri.configure(state="normal")
            self.cb_pot_pim.configure(state="normal")
            self.cb_pot_bmhftd.configure(state="normal")
        else:
            self.cb_pot_lj.configure(state=tk.DISABLED)
            self.cb_pot_ri.configure(state=tk.DISABLED)
            self.cb_pot_pim.configure(state=tk.DISABLED)
            self.cb_pot_bmhftd.configure(state=tk.DISABLED)
            self.control_values['lmorse'] = ".true."
            self.control_values['lvnlist'] = ".true."
            self.control_values['lreduced'] = ".true."

    #----------------------------------------------------------------------
    def _cb_pot_pim(self,event=None):
        if self.cv_pot['PIM'].get() == 0:
            self.cb_pot_morse.configure(state="normal")
            self.cb_pot_ri.configure(state="normal")
            self.cb_pot_lj.configure(state="normal")
            self.cb_pot_bmhftd.configure(state="normal")
        else:
            self.cb_pot_morse.configure(state=tk.DISABLED)
            self.cb_pot_ri.configure(state=tk.DISABLED)
            self.cb_pot_lj.configure(state=tk.DISABLED)
            self.cb_pot_bmhftd.configure(state=tk.DISABLED)
            self.field_values['lautoES'] = ".true."
            self.control_values['lcoulomb'] = ".true."
            self.control_values['lbmhftd'] = ".true."
            self.control_values['lbmhft'] = ".true."
            self.control_values['lreduced'] = ".true."

    #----------------------------------------------------------------------
    def _cb_pot_bmhftd(self,event=None):
        if self.cv_pot['BMHFTD'].get() == 0:
            self.cb_pot_morse.configure(state="normal")
            self.cb_pot_ri.configure(state="normal")
            self.cb_pot_pim.configure(state="normal")
            self.cb_pot_lj.configure(state="normal")
        else:
            self.cb_pot_morse.configure(state=tk.DISABLED)
            self.cb_pot_ri.configure(state=tk.DISABLED)
            self.cb_pot_pim.configure(state=tk.DISABLED)
            self.cb_pot_lj.configure(state=tk.DISABLED)

    #----------------------------------------------------------------------
    def _cb_integrator_nve(self,event=None):
        if self.cv_integrator['NVE'].get() == 0:
            self.cb_integrator_nvt.configure(state="normal")
            self.cb_integrator_npt.configure(state="normal")
        else:
            self.cb_integrator_nvt.configure(state=tk.DISABLED)
            self.cb_integrator_npt.configure(state=tk.DISABLED)
            self.md_values['integrator'] = "'nve-vv'"
    #----------------------------------------------------------------------
    def _cb_integrator_nvt(self,event=None):
        if self.cv_integrator['NVT'].get() == 0:
            self.cb_integrator_nve.configure(state="normal")
            self.cb_integrator_npt.configure(state="normal")
        else:
            self.cb_integrator_nve.configure(state=tk.DISABLED)
            self.cb_integrator_npt.configure(state=tk.DISABLED)
            self.md_values['integrator'] = "'nvt-nhcn'"
    #----------------------------------------------------------------------
    def _cb_integrator_npt(self,event=None):
        if self.cv_integrator['NPT'].get() == 0:
            self.cb_integrator_nve.configure(state="normal")
            self.cb_integrator_nvt.configure(state="normal")
        else:
            self.cb_integrator_nve.configure(state=tk.DISABLED)
            self.cb_integrator_nvt.configure(state=tk.DISABLED)
            self.md_values['integrator'] = "'npt-nhcnp'"


    #----------------------------------------------------------------------
    def _create_posff_file(self):
        # remove possible frame in the same region
        if self.plot_frame is not None :
            self.plot_frame.destroy()
        if self.control_md_frame is not None :    
            self.control_md_frame.destroy()

        self.config_frame = tk.Frame(self.parent, name='config_frame')
        self.config_frame.pack(side=tk.LEFT,fill=tk.BOTH)

        mainlabel = tk.Label(self.config_frame, text="Construct POSFF from scratch")
        mainlabel.grid(row=0,column=0,columnspan=5,sticky=tk.W+tk.E,pady=5,padx=10)

        self.config_frame.columnconfigure(0,weight=1)
        self.config_frame.columnconfigure(1,weight=1)
        self.config_frame.columnconfigure(2,weight=2)
        self.config_frame.columnconfigure(3,weight=2)
        self.config_frame.columnconfigure(4,weight=2)
        # nions
        lnions = tk.Label(self.config_frame, text="# ions")
        lnions.grid(row=1,column=3,sticky=tk.W+tk.E,pady=5,padx=5)
        self.nions=tk.IntVar(None)
        mnions=tk.Message(self.config_frame, textvariable=self.nions , justify=tk.LEFT, width=100)
        mnions.grid(row=1,column=4,columnspan=2,sticky=tk.W+tk.E,pady=5,padx=5)

        # atype
        latype = tk.Label(self.config_frame, text="ion type")
        latype.grid(row=2,column=3,sticky=tk.W+tk.E,pady=5,padx=5)
        self.atype=tk.StringVar(None)
        eatype=tk.Entry(self.config_frame,textvariable=self.atype)
        eatype.grid(row=2,column=4,columnspan=2,sticky=tk.W+tk.E,pady=5,padx=5)
        # itype
        litype = tk.Label(self.config_frame, text="nb per type")
        litype.grid(row=3,column=3,sticky=tk.W+tk.E,pady=5,padx=5)
        self.itype=tk.IntVar(None)
        eitype=tk.Entry(self.config_frame,textvariable=self.itype)
        eitype.grid(row=3,column=4,columnspan=2,sticky=tk.W+tk.E,pady=5,padx=5)

        # list box
        self.lb_type = MultiListbox(self.config_frame,(('type', 10),('nb per type', 10)))
        self.lb_type.grid(row=1,column=0,rowspan=4,columnspan=3,sticky=tk.W+tk.E,pady=5,padx=5)
        add_button = tk.Button(self.config_frame, text="Add", command = self._add_ion_to_list)
        add_button.grid(row=4,column=3,sticky=tk.W+tk.E,padx=10,pady=5)
        delete_button = tk.Button(self.config_frame, text="Remove",command=lambda lb_type=self.lb_type: self.lb_type.delete(tk.ACTIVE))
        delete_button.grid(row=4,column=4,sticky=tk.W+tk.E,pady=5,padx=10)


        self.cv_cubic = tk.IntVar()
        self.cv_ortho = tk.IntVar()
        self.cv_tric  = tk.IntVar()
        self.cv_cubic.set(1)
        self.cv_ortho.set(0)
        self.cv_tric.set(0)

        self.la = tk.Label(self.config_frame, text="a")
        self.la.grid(row=6,column=1,sticky=tk.W+tk.E,pady=5,padx=5)
        self.a1param=tk.StringVar(None)
        self.a1param.set(0.00000000)
        self.ea1=tk.Entry(self.config_frame,textvariable=self.a1param)
        self.ea1.grid(row=6,column=2,sticky=tk.W+tk.E,pady=5,padx=5)
        self.a2param=tk.StringVar(None)
        self.a2param.set(0.00000000)
        self.ea2=tk.Entry(self.config_frame,textvariable=self.a2param,state=tk.DISABLED)
        self.ea2.grid(row=6,column=3,sticky=tk.W+tk.E,pady=5,padx=5)
        self.a3param=tk.StringVar(None)
        self.a3param.set(0.00000000)
        self.ea3=tk.Entry(self.config_frame,textvariable=self.a3param,state=tk.DISABLED)
        self.ea3.grid(row=6,column=4,sticky=tk.W+tk.E,pady=5,padx=5)

        self.lb = tk.Label(self.config_frame, text="b",state=tk.DISABLED)
        self.lb.grid(row=7,column=1,sticky=tk.W+tk.E,pady=5,padx=5)
        self.b1param=tk.StringVar(None)
        self.b1param.set(0.00000000)
        self.eb1=tk.Entry(self.config_frame,textvariable=self.b1param,state=tk.DISABLED)
        self.eb1.grid(row=7,column=2,sticky=tk.W,pady=5,padx=5)
        self.b2param=tk.StringVar(None)
        self.b2param.set(0.00000000)
        self.eb2=tk.Entry(self.config_frame,textvariable=self.b2param,state=tk.DISABLED)
        self.eb2.grid(row=7,column=3,sticky=tk.W,pady=5,padx=5)
        self.b3param=tk.StringVar(None)
        self.b3param.set(0.00000000)
        self.eb3=tk.Entry(self.config_frame,textvariable=self.b3param,state=tk.DISABLED)
        self.eb3.grid(row=7,column=4,sticky=tk.W,pady=5,padx=5)

        self.lc = tk.Label(self.config_frame, text="c",state=tk.DISABLED)
        self.lc.grid(row=8,column=1,sticky=tk.W,pady=5,padx=5)
        self.c1param=tk.StringVar(None)
        self.c1param.set(0.00000000)
        self.ec1=tk.Entry(self.config_frame,textvariable=self.c1param,state=tk.DISABLED)
        self.ec1.grid(row=8,column=2,sticky=tk.W,pady=5,padx=5)
        self.c2param=tk.StringVar(None)
        self.c2param.set(0.00000000)
        self.ec2=tk.Entry(self.config_frame,textvariable=self.c2param,state=tk.DISABLED)
        self.ec2.grid(row=8,column=3,sticky=tk.W,pady=5,padx=5)
        self.c3param=tk.StringVar(None)
        self.c3param.set(0.00000000)
        self.ec3=tk.Entry(self.config_frame,textvariable=self.c3param,state=tk.DISABLED)
        self.ec3.grid(row=8,column=4,sticky=tk.W,pady=5,padx=5)

        self.cb_cubic = tk.Checkbutton(self.config_frame, text = "cubic", variable = self.cv_cubic, onvalue = 1, offvalue = 0 , command=self._cb_cubic)
        self.cb_cubic.grid(row=6,column=0,sticky=tk.W,pady=5)
        self.cb_ortho = tk.Checkbutton(self.config_frame, text = "ortorhombic", variable = self.cv_ortho, onvalue = 1, offvalue = 0 , command=self._cb_ortho,state=tk.DISABLED)
        self.cb_ortho.grid(row=7,column=0,sticky=tk.W,pady=5)
        self.cb_tric = tk.Checkbutton(self.config_frame, text = "triclinic", variable = self.cv_tric, onvalue = 1, offvalue = 0 , command=self._cb_tric,state=tk.DISABLED)
        self.cb_tric.grid(row=8,column=0,sticky=tk.W,pady=5)


        gene=tk.Button(self.config_frame, text="Generate random config", command = self._gene)
        gene.grid(row=9,column=0,columnspan=5,sticky=tk.W+tk.E,pady=5)

        filename='POSFF_'
    #----------------------------------------------------------------------
    def _gene(self):

        atype=[]
        itype=[]
        for l in self.lb_type.get(0,tk.END):
            atype.append(l[0]) 
            itype.append(l[1]) 

        if self.cv_cubic.get() == 1:
            u    =[ float(self.a1param.get()), 0.0         , 0.0          ]
            v    =[ 0.0         , float(self.a1param.get()), 0.0          ]
            w    =[ 0.0         , 0.0         , float(self.a1param.get()) ]
        if self.cv_ortho.get() == 1:
            u    =[ float(self.a1param.get()), 0.0         , 0.0          ]
            v    =[ 0.0         , float(self.b1param.get()), 0.0          ]
            w    =[ 0.0         , 0.0         , float(self.c1param.get()) ]
        elif self.cv_tric.get() == 1:  
            u    =[ float(self.a1param.get()) , float(self.a2param.get()) ,float(self.a3param.get()) ]
            v    =[  float(self.b1param.get()) , float(self.b2param.get()) ,float(self.b3param.get()) ]
            w    =[  float(self.c1param.get()) , float(self.c2param.get()) ,float(self.c3param.get()) ]
       
        ions=[]
        conf = Config(ions=ions,u=u,v=v,w=w,system='randomPY',natm=self.nions.get(),ntype=len(atype),types=atype, natmpertype = itype, coord_format='Direct')
        
        for ia in range(conf.natm):
                conf.ions.append (  Ion ( index_ion=ia ) )

        conf.typeinfo_init()

        # random structure in the box
        for ia in range(conf.natm):
            x = random.random()
            y = random.random()
            z = random.random()
            conf.move(ia,pos=[x,y,z])

        conf.write_POSFF('POSFF.randomPY')
        file = open('POSFF.randomPY', 'rb')
        if file != None:
            content = file.read()
            self.txt_posff.delete('1.0','end')
            self.txt_posff.insert('1.0',content)
            file.close()

        self.nb.select(1) # show tab 1 == posff

    #----------------------------------------------------------------------
    def _cb_cubic(self,event=None):
        if self.cv_cubic.get() == 0:
            self.cb_ortho.configure(state="normal")
            self.cb_tric.configure(state="normal")
            self.la.configure(state=tk.DISABLED)
            self.lb.configure(state=tk.DISABLED)
            self.lc.configure(state=tk.DISABLED)
            self.ea1.configure(state=tk.DISABLED)
            self.ea2.configure(state=tk.DISABLED)
            self.ea3.configure(state=tk.DISABLED)
            self.eb1.configure(state=tk.DISABLED)
            self.eb2.configure(state=tk.DISABLED)
            self.eb3.configure(state=tk.DISABLED)
            self.ec1.configure(state=tk.DISABLED)
            self.ec2.configure(state=tk.DISABLED)
            self.ec3.configure(state=tk.DISABLED)
        else:
            self.cb_ortho.configure(state=tk.DISABLED)
            self.cb_tric.configure(state=tk.DISABLED)
            self.la.configure(state="normal")
            self.lb.configure(state=tk.DISABLED)
            self.lc.configure(state=tk.DISABLED)
            self.ea1.configure(state="normal")
            self.ea2.configure(state=tk.DISABLED)
            self.ea3.configure(state=tk.DISABLED)
            self.eb1.configure(state=tk.DISABLED)
            self.eb2.configure(state=tk.DISABLED)
            self.eb3.configure(state=tk.DISABLED)
            self.ec1.configure(state=tk.DISABLED)
            self.ec2.configure(state=tk.DISABLED)
            self.ec3.configure(state=tk.DISABLED)


    #----------------------------------------------------------------------
    def _cb_ortho(self,event=None):
        if self.cv_ortho.get() == 0:
            self.cb_cubic.configure(state="normal")
            self.cb_tric.configure(state="normal")
            self.la.configure(state=tk.DISABLED)
            self.lb.configure(state=tk.DISABLED)
            self.lc.configure(state=tk.DISABLED)
            self.ea1.configure(state=tk.DISABLED)
            self.ea2.configure(state=tk.DISABLED)
            self.ea3.configure(state=tk.DISABLED)
            self.eb1.configure(state=tk.DISABLED)
            self.eb2.configure(state=tk.DISABLED)
            self.eb3.configure(state=tk.DISABLED)
            self.ec1.configure(state=tk.DISABLED)
            self.ec2.configure(state=tk.DISABLED)
            self.ec3.configure(state=tk.DISABLED)
        else:
            self.cb_cubic.configure(state=tk.DISABLED)
            self.cb_tric.configure(state=tk.DISABLED)
            self.la.configure(state="normal")
            self.lb.configure(state="normal")
            self.lc.configure(state="normal")
            self.ea1.configure(state="normal")
            self.ea2.configure(state=tk.DISABLED)
            self.ea3.configure(state=tk.DISABLED)
            self.eb1.configure(state="normal")
            self.eb2.configure(state=tk.DISABLED)
            self.eb3.configure(state=tk.DISABLED)
            self.ec1.configure(state="normal")
            self.ec2.configure(state=tk.DISABLED)
            self.ec3.configure(state=tk.DISABLED)
    #----------------------------------------------------------------------
    def _cb_tric(self,event=None):
        if self.cv_tric.get() == 0:
            self.cb_ortho.configure(state="normal")
            self.cb_cubic.configure(state="normal")
            self.la.configure(state=tk.DISABLED)
            self.lb.configure(state=tk.DISABLED)
            self.lc.configure(state=tk.DISABLED)
            self.ea1.configure(state=tk.DISABLED)
            self.ea2.configure(state=tk.DISABLED)
            self.ea3.configure(state=tk.DISABLED)
            self.eb1.configure(state=tk.DISABLED)
            self.eb2.configure(state=tk.DISABLED)
            self.eb3.configure(state=tk.DISABLED)
            self.ec1.configure(state=tk.DISABLED)
            self.ec2.configure(state=tk.DISABLED)
            self.ec3.configure(state=tk.DISABLED)
        else:
            self.cb_cubic.configure(state=tk.DISABLED)
            self.cb_ortho.configure(state=tk.DISABLED)
            self.la.configure(state="normal")
            self.lb.configure(state="normal")
            self.lc.configure(state="normal")
            self.ea1.configure(state="normal")
            self.ea2.configure(state="normal")
            self.ea3.configure(state="normal")
            self.eb1.configure(state="normal")
            self.eb2.configure(state="normal")
            self.eb3.configure(state="normal")
            self.ec1.configure(state="normal")
            self.ec2.configure(state="normal")
            self.ec3.configure(state="normal")

    #----------------------------------------------------------------------
    def _add_ion_to_list(self):
        if self.atype.get() != "" and self.itype.get() != 0:
            print "adding type",self.atype.get(),self.itype.get()
            self.lb_type.insert (tk.END, (self.atype.get(),self.itype.get()))
            self.nions.set(0)
            for l in self.lb_type.get(0,tk.END):
                self.nions.set(self.nions.get()+l[1])
            print "current number of total ions ",self.nions.get()

    #----------------------------------------------------------------------
    def _apropos(self):
        apropos = tk.Toplevel()
        apropos.resizable(width=tk.FALSE, height=tk.FALSE)
        apropos.title("About MDff")
        about_message="""
        MDFF parallel Molecular Dynamics ... For Fun
        Copyright (C) 2011  F. Vasconcelos
        
        This program is free software; you can redistribute it and/or
        modify it under the terms of the GNU General Public License
        as published by the Free Software Foundation; either version 2
        of the License, or (at your option) any later version.
        
        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.
        
        You should have received a copy of the GNU General Public License
        along with this program; if not, WRITE to the Free Software
        Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
        02110-1301, USA.
        """
        msg = tk.Message(apropos, text=about_message , justify=tk.LEFT, relief=tk.RAISED,width=2000)
        msg.pack()

    #----------------------------------------------------------------------
    def _doc(self):
        pdf_file_path = self.doc_path + '/' + 'mdff.pdf'

        sys_pltfrm = platform.system()
        print sys_pltfrm
        cmd = { 'Linux'  :'xpdf '}

        try:
            subprocess.call([cmd[sys_pltfrm] + pdf_file_path], shell=True)
        except OSError as err:
            sys.exit('Error opening pdf:\n\t{}'.format(err))


    #----------------------------------------------------------------------
    def _save_as(self):
        if self.nb.tab(self.nb.select(), "text") == "control.F" :
            name = tkFileDialog.asksaveasfile(mode='w',filetypes=[('control files','.F'),('all files','.*')])
            text2save=str(self.txt_control.get(0.0,tk.END))
            name.write(text2save)
            name.close
        elif self.nb.tab(self.nb.select(), "text") == "POSFF":
            name = tkFileDialog.asksaveasfile(mode='w',filetypes=[('POSFF','POSFF* CONTFF*'),('all files','.*')])
            text2save=str(self.txt_posff.get(0.0,tk.END))
            name.write(text2save)
            name.close
        elif self.nb.tab(self.nb.select(), "text") == "log":
            name = tkFileDialog.asksaveasfile(mode='w',filetypes=[('all files','.*')])
            text2save=str(self.txt_stdout.get(0.0,tk.END))
            name.write(text2save)
            name.close

    #----------------------------------------------------------------------
    def _select_control_file(self):
        self.control_path.set(tkFileDialog.askopenfilename(title="Select MDff control file",filetypes=[('control files','.F'),('all files','.*')]))
        file = open(self.control_path.get(), 'a+')
        if file != None:
            content = file.read()
            self.txt_control.delete('1.0','end')
            self.txt_control.insert('1.0',content)
            file.close()
        self.txt_control.update()
        self.nb.select(0) # show control tab

    #----------------------------------------------------------------------
    def _select_posff_file(self):
        self.posff_path.set(tkFileDialog.askopenfilename(title="Select MDff POSFF file",filetypes=[('POSFF files','POSFF* CONTFF*'),('all files','*')]))
        file = open(self.posff_path.get(), 'a+')
        if file != None:
            content = file.read()
            self.txt_posff.delete('1.0','end')
            self.txt_posff.insert('1.0',content)
            file.close()
        self.nb.select(1) # show posff tab
    #----------------------------------------------------------------------
    def _select_exe(self):
        self.exe_path.set(tkFileDialog.askopenfilename(title="Select MDff executable",filetypes=[('exe files','.x'),('all files','.*')]))
    #----------------------------------------------------------------------
    def start_exe(self,event):
        self._start_exe()
    #----------------------------------------------------------------------
    def _start_exe(self):

        os.system("cp "+self.posff_path.get()+" POSFF")
        self.nb.select(2) # show log tab
        self.proc = subprocess.Popen('mpirun '+self.exe_path.get()+' '+self.control_path.get(), stdout=subprocess.PIPE ,shell=True, preexec_fn=os.setsid)
        while True:
            line = self.proc.stdout.readline()
            if line != '':
                l = line.rstrip()
                self.txt_stdout.insert(tk.END, l+"\n")
                self.txt_stdout.see(tk.END) 
                self.txt_stdout.update()

            else:
                break
        

    #----------------------------------------------------------------------
    def _kill_exe(self,event):

        os.killpg(os.getpgid(self.proc.pid), signal.SIGTERM)

    #----------------------------------------------------------------------
    def _quit(self,event):
        sys.exit(0)
    #----------------------------------------------------------------------
    def _save(self,event):
        if self.nb.tab(self.nb.select(), "text") == "control.F" :
            if os.path.isfile(self.control_path.get()) :
                file_ = open(self.control_path.get(),'w')
                text2save=str(self.txt_control.get(0.0,tk.END))
                file_.write(text2save)
                file_.close
            else:
                self._save_as()
        elif self.nb.tab(self.nb.select(), "text") == "POSFF":
            if os.path.isfile(self.posff_path.get()) :
                file_ = open(self.posff_path,'w')
                text2save=str(self.txt_posff.get(0.0,tk.END))
                file_.write(text2save)
                file_.close
            else:
                self._save_as()
        elif self.nb.tab(self.nb.select(), "text") == "log":
            if os.path.isfile(self.stdout_path.get()) :
                file_ = open(self.stdout_path.get(),'w')
                text2save=str(self.txt_stdout.get(0.0,tk.END))
                file_.write(text2save)
                file_.close
            else:
                self._save_as()
    #----------------------------------------------------------------------
    def _run_stat(self):
        self.tstat = tk.Toplevel()
        self.tstat.title("How many points ?")
        lstat = tk.Label(self.tstat,text="#points ? (all)")
        lstat.grid(row=0,column=0)
        self.last_points=tk.StringVar(None)
        estat = tk.Entry(self.tstat,textvariable=self.last_points)
        estat.grid(row=1,column=0)
        #estat.bind("<Return>",self._kill_window())
        self.tstat.wait_window()
        self._run_poszi()

    #----------------------------------------------------------------------
    def _kill_window(self):
        print "inside _kill_window"
        self.tstat.quit()

    #----------------------------------------------------------------------
    def _run_poszi(self):

        filename='OSZIFF'
        if not os.path.isfile('OSZIFF') :
            tkMessageBox.showinfo("File missing", "OSZIFF file is missing for analysis")
            return
        if format_filename(filename):
            self.name_quant,self.alldata = read_OSZIFF(filename)
        else:
            raise ValueError('Format of input file doesnt match OSZIFF')
     
        if self.last_points.get() == 'all' or self.last_points.get() == '' :
            lp = len(self.alldata[0])
        else:
            lp = int( self.last_points.get() )

        print lp,self.last_points.get()
    
        self.average,self.ainfo=averaging(self.name_quant,self.alldata,lp)

        self._create_plot() 
        self._average_info()
    
    #----------------------------------------------------------------------
    def _new_plot(self,lp):

        self.f=Figure()
        self.a = self.f.add_subplot(111)
        q=[2,12]
        title='Total Energy MD'
        t = self.alldata[1][-lp:]
        s = self.alldata[q[0]][-lp:]
        x1= [self.average[q[0]]]*len(t)
        self.a.plot(t, s,'b')
        self.a.plot(t, x1, '--' , color='b')
        s = self.alldata[q[1]][-lp:]
        x1= [self.average[q[1]]]*len(t)
        self.a.plot(t, s , 'g')
        self.a.plot(t, x1, '--', color='g')
        self.a.set_xlabel(self.name_quant[1]+units[self.name_quant[1]])
        self.a.set_ylabel(self.name_quant[q[0]]+units[self.name_quant[q[0]]])
        self.a.set_title(title)

    #----------------------------------------------------------------------
    def _update_plot(self,delta,lp):
        self.a.clear()
        self.k+=delta
        if self.k > 3 :
            self.k = 3
            return
        if self.k < 0 :
            self.k = 0
            return
        self.a.clear()
        titles=['Total Energy MD','Potential Energy MD','Temperature MD','Volume MD']
        qs=[2,4,7,11]
        q=qs[self.k]
        title=titles[self.k]
        t = self.alldata[1][-lp:]
        s = self.alldata[q][-lp:]
        x1= [self.average[q]]*len(t)
        self.a.plot(t, s,'b')
        self.a.plot(t, x1, '--' , color='b')
        if self.k == 0 :
            s = self.alldata[12][-lp:]
            x1= [self.average[12]]*len(t)
            self.a.plot(t, s , 'g')
        self.a.set_xlabel(self.name_quant[1]+units[self.name_quant[1]])
        self.a.set_ylabel(self.name_quant[q]+units[self.name_quant[q]])
        self.a.set_title(title)
        self.canvas.draw()

    #----------------------------------------------------------------------
    def _create_plot(self):
        if self.config_frame is not None :
            self.config_frame.destroy()
        if self.control_md_frame is not None :    
            self.control_md_frame.destroy()
        self.plot_frame = tk.Frame(self.parent, name='plot_frame')
        self.plot_frame.pack(side=tk.LEFT,fill=tk.Y)


        if self.last_points.get() == 'all' or self.last_points.get() == '' :
            lp = len(self.alldata[0])
        else:
            lp = int( self.last_points.get() )
        self._new_plot(lp)
        self.canvas = FigureCanvasTkAgg(self.f, self.plot_frame)
        self.canvas.draw()
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.X)

        toolbar = NavigationToolbar2TkAgg(self.canvas, self.plot_frame)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP,fill=tk.X)

        self.k=0
        previous_button = tk.Button(self.plot_frame, text="Previous Plot", command = lambda : self._update_plot(-1,lp),width=10)
        next_button = tk.Button(self.plot_frame, text="Next Plot", command = lambda : self._update_plot(1,lp),width=10)
        previous_button.pack(side=tk.LEFT)
        next_button.pack(side=tk.LEFT)
    #----------------------------------------------------------------------
    def _average_info(self):
        ainfo = tk.Toplevel()
        ainfo.title("Run statistics")
        ainfo.geometry('%dx%d+%d+%d' % (700, 320, 0, 0) )
        ainfo.update()
        ainfo_message='\n'
        ainfo_message=str(len(self.alldata[0]))+" points in input file\n"
        ainfo_message+="averaging on last "+str(self.last_points.get())+" points\n"
        ainfo_message+='\n'
        for l in self.ainfo:
            ll='{0:<14} {1:<10} {2:15.8e} {3:^10} {4:15.8e} \n'.format("<"+l[0]+">","=",l[1],"std.",l[2])
            ll.ljust(70)
            ainfo_message+=ll
        font = tkFont.Font(family="Courier", size=11,weight='bold')
        msg = tk.Message(ainfo, text=ainfo_message , anchor='center',width=1000,font=font)
        msg.pack(expand=True, fill='x',padx=10,pady=20)
        msg.bind("<Configure>", lambda e: msg.configure(width=e.width-10))

    #----------------------------------------------------------------------
    def _nothing(self):
        pass

if __name__ == "__main__":

    root = tk.Tk()
    app  = MApp(root)
    root.mainloop()

