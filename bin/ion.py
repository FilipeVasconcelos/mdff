# =====================================================================================
import numpy as np
import os
import sys
import time
import random

from lattice import Lattice
from vesta import vesta_definitions
# =====================================================================================

# =====================================================================================
class Ion(object) :
    """ 
    Definition of an ion:
    ------
    In this implementation, we defined each ion as an object 
    and Config as an list of ions.
    Note :
    Their is several operations (methods) on this object. 
    I try to follow the golden rule : "only change attributs through methods"
    """
    
    # ==========================================================================
    def __init__(self,index_type=0,type='',index_ion=0,\
                 pos=[0.,0.,0.],vel=[0.,0.,0.],force=[0.,0.,0.],\
                 charge=0.0,dipole=[0.,0.,0.],quadmom=0.0,allowedmove=[ 'T' , 'T',  'T' ],\
                 diso=0.0 , rgb =None , radius = None):

        self.index_ion     = index_ion
        self.type          = type 
        self.index_type    = index_type 
        self.r             = pos 
        self.v             = vel
        self.f             = force
        self.charge        = charge
        self.dipole        = dipole 
        self.quadia        = quadmom
        self.label         = type+str(index_ion+1)  
        self.allowedmove   = allowedmove
    
        # nmr related :
        self.diso          = diso
        self.labelnmr      = str(round(float(self.diso),2))+'ppm' 

        # radius and rgb data for gen_vesta
        if radius == None and type != '' :
            self.radius , dumb_rgb = vesta_definitions(type)
        else :
            self.radius  = radius
        if rgb == None and type != '' :
            dumb_radius , self.rgb = vesta_definitions(type) 
        else :
            self.rgb     = rgb

    # =====================================================================================
    def zero_forces(self):
        self.f = [ 0., 0. , 0. ]
    # =====================================================================================
    def zero_velocities(self):
        self.v = [ 0., 0. , 0. ]
    # =====================================================================================
    def change_typeinfo(self, type , index_type ):
        self.type        = type
        self.label       = type+str(self.index_ion)
        self.index_type  = index_type
        self.radius, self.rgb     = vesta_definitions(type)
        return self.type, self.label, self.index_type , self.radius , self.rgb
    # =====================================================================================
    def change_nmrinfo(self, diso ):
        self.diso          = diso
        self.labelnmr      = str(round(float(self.diso),2))+'ppm'
        return self.diso, self.labelnmr
    # =====================================================================================
    def get_coord(self):
        return self.r
    # =====================================================================================
    def get_velocity(self):
        return self.v
    # =====================================================================================
    def get_forces(self):
        return self.f
    # =====================================================================================
    def move (self,pos):
        return [ self.r[0] + pos[0], self.r[1] + pos[1], self.r[2] + pos[2] ]  
    # =====================================================================================
    def add_to_velocity(self,vel):
        return [ self.v[0] + vel[0], self.v[1] + vel[1], self.v[2] + vel[2] ]  
    # =====================================================================================
    def add_to_forces(self,f):
        return [ self.f[0] + f[0], self.f[1] + f[1], self.f[2] + f[2] ]  
    # =====================================================================================
    def print_ion(self):
        pass


if __name__ == '__main__':

    u = [ 10.0,  0.,  0. ]
    v = [  0.0, 10.,  0. ]
    w = [  0.0,  0., 10. ]
    
    ions=[]
    conf_test=Config(ions=ions,u=u,v=v,w=w)

    print( conf_test.print_direct_basis(None) )

