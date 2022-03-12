# =====================================================================================
import numpy as np
import os
import sys
import time
import random

from ion import Ion
from lattice import Lattice
from constants import float1,float3,posff,poscar,config
from vesta import vesta_definitions
# =====================================================================================
class Config(Lattice):
    """
    Configuration Class :
    -------
    it has Lattice class as parent class and list of ions. 
    Note :
    (i)  It seems that I could have consider the class Ion as 
         an subclass of Config ... but it works fine it this way. 
    (ii) Several methods of this class ( change_typeinfo, change_nmrinfo) are "drivers" methods 
         to Ion methods with the same name.
         Some are just loops of Ion methods for all ions ( method_to_all ) 
    """

    # =====================================================================================
    def __init__( self, ions =[], system='', nion=0, ntype=0, types=[],\
                  natmpertype = [], coord_format='C', rcut_coeff = 1.0,\
                  pbc = True , u=[1.0,0.0,0.0], v=[0.0,1.0,0.0], w=[0.0,0.0,1.0]):
        self.ions               = ions
        self.system             = system
        self.nion               = nion 
        self.ntype              = ntype
        self.types              = types
        self.natmpertype        = natmpertype
        self.coord_format       = coord_format
        self.dyn                = False 
        self.cartesian_allowed  = ['C', 'Cartesian' ]
        self.direct_allowed     = ['D', 'Direct' ]
        self.rcut_coeff         = rcut_coeff
        self.pbc	        = pbc
        super(Config, self).__init__(u=u,v=v,w=w)  # construct __init__ method of the superclass (Lattice)

    # =====================================================================================
    def change_typeinfo(self, index_ion, type, index_type ):
        self.ions[index_ion].type, 
        self.ions[index_ion].label, 
        self.ions[index_ion].index_type, 
        self.ions[index_ion].radius, 
        self.ions[index_ion].rgb = self.ions[index_ion].change_typeinfo( type , index_type )

    # =====================================================================================
    def change_nmrinfo(self, index_ion, diso ):
        self.ions[index_ion].diso  , self.ions[index_ion].labelnmr = self.ions[index_ion].change_nmrinfo( diso )
	
    # =====================================================================================
    def print_sample(self):
        pass
    # =====================================================================================
    def zero_forces_to_all(self):
        for item in self.ions:
            item.zero_forces()
    # =====================================================================================
    def zero_velocities_to_all(self):
        for item in self.ions:
            item.zero_velocities()
    # =====================================================================================
    def move (self,ia,pos):
        self.ions[ia].r = self.ions[ia].move ( pos )
    # =====================================================================================
    def add_to_velocity(self,ia,vel):
        self.ions[ia].v = self.ions[ia].add_to_velocity( vel )
    # =====================================================================================
    def add_to_forces(self,ia,force):
        self.ions[ia].f = self.ions[ia].add_to_forces( force )

    # =====================================================================================
    def calc_temp(self):
        ekin = 0.0
        for ia in range ( self.nion ) :
            vx = self.ions[ia].v[0]
            vy = self.ions[ia].v[1]
            vz = self.ions[ia].v[2]
            vx2 = vx * vx 
            vy2 = vy * vy 
            vz2 = vz * vz 
            ekin +=  vx2 + vy2 + vz2 
            ekin = ekin * 0.5
            T    = ( 2. / 3. ) * ekin
            T    /= float ( self.nion ) 
            return T , ekin

    # =====================================================================================
    def dist_2atoms(self,ia,ja):
        rs = self.ions[ia].get_coord()
        ro = self.ions[ja].get_coord()
        dx = rs[0] - ro[0]
        dy = rs[1] - ro[1]
        dz = rs[2] - ro[2]
        sx = dx - round ( dx )
        sy = dy - round ( dy )
        sz = dz - round ( dz )
        dx = sx * self.direct[0][0] + sy * self.direct[0][1] + sz * self.direct[0][2]
        dy = sx * self.direct[1][0] + sy * self.direct[1][1] + sz * self.direct[1][2]
        dz = sx * self.direct[2][0] + sy * self.direct[2][1] + sz * self.direct[2][2]
        dist_square = dx * dx + dy * dy + dz * dz
        self.vector_ij = []
        self.vector_ij = [ dx , dy , dz ]
        return dist_square
    # =====================================================================================
    def dist_2atoms_nopbc(self,ia,ja):
       	rs = self.ions[ia].get_coord()
        ro = self.ions[ja].get_coord()
        sx = rs[0] - ro[0]
        sy = rs[1] - ro[1]
        sz = rs[2] - ro[2]
        dx = sx * self.direct[0][0] + sy * self.direct[0][1] + sz * self.direct[0][2]
        dy = sx * self.direct[1][0] + sy * self.direct[1][1] + sz * self.direct[1][2]
        dz = sx * self.direct[2][0] + sy * self.direct[2][1] + sz * self.direct[2][2]
        dist_square = dx * dx + dy * dy + dz * dz
        self.vector_ij = []
        self.vector_ij = [ dx , dy , dz ]
        return dist_square

    # =====================================================================================
    def periodicbc(self):

        for ia in range(self.nion):
            self.ions[ia][0] = self.ions[ia].r[0] - round ( self.ions[ia][0] ) 			
            self.ions[ia][1] = self.ions[ia].r[1] - round ( self.ions[ia][1] ) 			
            self.ions[ia][2] = self.ions[ia].r[2] - round ( self.ions[ia][2] ) 			

    # =====================================================================================
    def typeinfo_init(self):    
        cc = 0	
        for it in range(self.ntype):
            ccs = cc
            cc = cc + self.natmpertype[it]
            for ia in range ( ccs , cc  ):
                self.change_typeinfo ( ia ,  type =self.types[it] , index_type = it )


    # =====================================================================================
    def kardir ( self ):
        """
        transform a set of vectors from cartesian coordinates to
        ) direct lattice      (BASIS must be equal to B reciprocal lattice)
        ) reciprocal lattice  (BASIS must be equal to A direct lattice)
        adapted from vasp
        note here only cartesian coordinates --> direct lattice
        """
        if self.coord_format not in self.cartesian_allowed : 
            print( "coordinates already in Direct",self.coord_format,self.cartesian_allowed)
            raise ValueError('in Config.kardir ')
        for ia in range( len ( self.ions ) ) :
            v1 = self.ions[ia].r[0] * self.reciprocal[0][0] + \
                 self.ions[ia].r[1] * self.reciprocal[0][1] + \
                 self.ions[ia].r[2] * self.reciprocal[0][2]
            v2 = self.ions[ia].r[0] * self.reciprocal[1][0] + \
                 self.ions[ia].r[1] * self.reciprocal[1][1] + \
                 self.ions[ia].r[2] * self.reciprocal[1][2]
            v3 = self.ions[ia].r[0] * self.reciprocal[2][0] + \
                 self.ions[ia].r[1] * self.reciprocal[2][1] + \
                 self.ions[ia].r[2] * self.reciprocal[2][2]
            self.ions[ia].r[0]=v1
            self.ions[ia].r[1]=v2
            self.ions[ia].r[2]=v3
            self.coord_format = 'Direct'

    # =====================================================================================
    def dirkar ( self ):
        """
        transform a set of vectors from
        ) direct lattice      (BASIS must be equal to A direct lattice)
        ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
         to cartesian coordinates
         adapted from VASP
        note here only direct lattice --> cartesian coordinates
        """
        if self.coord_format not in self.direct_allowed : 
            print( "coordinates already in Cartesian")
            raise ValueError('in Config.dirkar ')
        for ia in range( len ( self.ions ) ) :
            v1 = self.ions[ia].r[0] * self.direct[0][0] + \
                 self.ions[ia].r[1] * self.direct[1][0] + \
                 self.ions[ia].r[2] * self.direct[2][0]
            v2 = self.ions[ia].r[0] * self.direct[0][1] + \
                 self.ions[ia].r[1] * self.direct[1][1] + \
                 self.ions[ia].r[2] * self.direct[2][1]
            v3 = self.ions[ia].r[0] * self.direct[0][2] + \
                 self.ions[ia].r[1] * self.direct[1][2] + \
                 self.ions[ia].r[2] * self.direct[2][2]
            self.ions[ia].r[0]=v1
            self.ions[ia].r[1]=v2
            self.ions[ia].r[2]=v3
            self.coord_format = 'Cartesian'

    # =====================================================================================
    def change_rcut_coeff(self,rcut) :
        self.rcut_coeff += rcut 
        if self.rcut_coeff > 1. :
            self.rcut_coeff = 1.    
        dexp = self.exp_dist()
        return dexp

    # =====================================================================================
    def read_POSFF(self,filename):
        try:
            f = open(filename)
        except IOError:
            sys.exit()
        else:
            lines       = f.readlines()
        
        self.nion   = int(lines[0])
        self.system = lines[1].strip()
        u           = [ float(item) for item in lines[2].split() ]
        v           = [ float(item) for item in lines[3].split() ]
        w           = [ float(item) for item in lines[4].split() ]
        assert len(u)==3 , "line 3 of POSFF should have only 3 elements"
        assert len(v)==3 , "line 4 of POSFF should have only 3 elements"
        assert len(w)==3 , "line 5 of POSFF should have only 3 elements"

        self.ntype        = int(lines[5])
        self.types        = [ item for item in lines[6].split() ]
        self.natmpertype  = [ int(item) for item in lines[7].split() ]
        self.coord_format = lines[8].strip()

        for i in range(9, 9+self.nion):
            line_list = [ item for item in lines[i].split()]
            self.ions.append (Ion(index_ion=i-8, \
                              type=line_list[0], \
                              pos = [ float(line_list[1]), float(line_list[2]) , float(line_list[3]) ]))

        self.typeinfo_init()		
        self.lattice_reinit(u,v,w,1.0)
        self.volumeANDoppositefaces(self.direct)
        self.lattice_parameters (self.direct)
        self.reciprocal_basis (self.direct)

        f.close()
    
    # =====================================================================================
    def write_XYZ(self,filename):
        f = open(filename,'w+')
        f.write(str(self.nion)+'\n')
        f.write(self.system+'\n')
        for ia in range(self.nion):
            f.write(posff.format(str(self.ions[ia].type),\
                                        self.ions[ia].r[0],\
                                        self.ions[ia].r[1],\
                                        self.ions[ia].r[2]))
        f.close()

    # =====================================================================================
    def write_CONFIG(self,filename):
        f = open(filename,'w+')
        f.write(self.system+'\n')
        f.write("0  1"+str(self.nion)+" 0.0 \n")
        self.print_direct_basis(f) # standard output filename = None
        for ia in range(self.nion):
            f.write(self.ions[ia].type+"         "+str(ia)+"\n")
            f.write(config.format(self.ions[ia].r[0],self.ions[ia].r[1],self.ions[ia].r[2]))
        f.close()
    # =====================================================================================
    def write_POSFF(self,filename):

        f = open(filename,'w+')
        f.write(str(self.nion)+'\n')
        f.write(self.system+'\n')
        self.print_direct_basis(f) # standard output filename = None
        f.write(str(self.ntype)+'\n')
        for k in range(self.ntype-1):
            f.write(str(self.types[k])+" ")
        f.write(str(self.types[-1])+"\n")
        for k in range(self.ntype-1):
            f.write(str(self.natmpertype[k])+" ")
        f.write(str(self.natmpertype[-1])+"\n")
        f.write(self.coord_format+'\n')
        for ia in range(self.nion):
            f.write(posff.format(str(self.ions[ia].type),\
                                        self.ions[ia].r[0],\
                                        self.ions[ia].r[1],\
                                        self.ions[ia].r[2]))
        f.close()
    # =====================================================================================
    def write_POSCAR(self,filename):

        f = open(filename,'w+')
        f.write(self.system+'\n')
        f.write(str(self.coeff)+'\n')
        self.print_direct_basis(f) # standard output filename = None
        for k in range(self.ntype-1):
            f.write(str(self.types[k])+" ")
        f.write(str(self.types[-1])+"\n")
        for k in range(self.ntype-1):
            f.write(str(self.natmpertype[k])+" ")
        f.write(str(self.natmpertype[-1])+"\n")
        f.write(self.coord_format+'\n')
        for ia in range(self.nion):
            f.write(poscar.format(self.ions[ia].r[0],\
                                  self.ions[ia].r[1],\
                                  self.ions[ia].r[2],str(self.ions[ia].type)))
        f.close()
    # =====================================================================================
    def write_TRAJFF(self,filename,option):
    
        pass

# ============================
#        testing unit  :
# ============================
if __name__ == '__main__':

    u = [ 10.0,  0.,  0. ]
    v = [  0.0, 10.,  0. ]
    w = [  0.0,  0., 10. ]
    
    ions=[]
    conf_test=Config(ions=ions,u=u,v=v,w=w)
    conf_test.print_direct_basis(None)
    conf_test.read_POSFF('POSFF')
    conf_test.write_POSFF('POSFF.test')

