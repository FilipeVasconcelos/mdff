# =====================================================================================
import numpy as np
import os
import sys
import time
import random

from lattice import Lattice
from constants import float_format
from vesta import vesta_definitions
# =====================================================================================

# =====================================================================================
class Ion(object) :
    """ 
    Definition of an ion:
    ------
    In this implementation, we defined each ion as an object and Config as an list of ions.
    Note :
    Their is several operations (methods) on this object. 
    I try to follow the golden rule : "only change attributs through methods"
    """
    
    # =====================================================================================
    def __init__(self,index_type=0,type='',index_ion=0,pos=[0.,0.,0.],vel=[0.,0.,0.],force=[0.,0.,0.],charge=0.0,dipole=[0.,0.,0.],quadmom=0.0,allowedmove=[ 'T' , 'T',  'T' ] , diso=0.0 , rgb =None , radius = None):

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
	print '%4s' % str(self.ions[i].label) , float_format % self.ions[i].qia , form % ( self.ions[i].r[0], self.ions[i].v[0], self.ions[i].f[0] )



# ===============================================================================================================================================================
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
    def __init__( self, ions =[], system='', natm=0, ntype=0, types=[], natmpertype = [], coord_format='C', rcut_coeff = 1.0 ,pbc = True , u=[1.0,0.0,0.0], v=[0.0,1.0,0.0], w=[0.0,0.0,1.0]):
        self.ions               = ions
	self.system             = system
	self.natm               = natm
	self.ntype              = ntype
	self.types              = types
	self.natmpertype        = natmpertype
	self.coord_format       = coord_format
	self.dyn                = True
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

	form  = 4*float_format
        print 80*'-'
        print "                             sample of config"
        print 80*'-'
        print "%4s%16s%16s%16s%16s" % ("site","rx","ry","rz","diso")
	if self.natm < 13 :
	    for i in range(self.natm):
	        print '%4s' % str(self.ions[i].label) , form % ( self.ions[i].r[0], self.ions[i].r[1], self.ions[i].r[2] , self.ions[i].diso )
	else:
            for i in range(4):	
	        print '%4s' % str(self.ions[i].label) , form % ( self.ions[i].r[0], self.ions[i].r[1], self.ions[i].r[2] , self.ions[i].diso )
	    for i in range(self.natm-4,self.natm):
	        print '%4s' % str(self.ions[i].label) , form % ( self.ions[i].r[0], self.ions[i].r[1], self.ions[i].r[2] , self.ions[i].diso )


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
	for ia in range ( self.natm ) :
	    vx = self.ions[ia].v[0]
	    vy = self.ions[ia].v[1]
	    vz = self.ions[ia].v[2]
	    vx2 = vx * vx 
	    vy2 = vy * vy 
	    vz2 = vz * vz 
	    ekin +=  vx2 + vy2 + vz2 
  
	    ekin = ekin * 0.5
	    T    = ( 2. / 3. ) * ekin
	    T    /= float ( self.natm ) 
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

        for ia in range(self.natm):
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
	    print "coordinates already in Direct",self.coord_format,self.cartesian_allowed
	    raise ValueError('in Config.kardir ')
	for ia in range( len ( self.ions ) ) :
	    v1 = self.ions[ia].r[0] * self.reciprocal[0][0] + self.ions[ia].r[1] * self.reciprocal[0][1] + self.ions[ia].r[2] * self.reciprocal[0][2]
	    v2 = self.ions[ia].r[0] * self.reciprocal[1][0] + self.ions[ia].r[1] * self.reciprocal[1][1] + self.ions[ia].r[2] * self.reciprocal[1][2]
	    v3 = self.ions[ia].r[0] * self.reciprocal[2][0] + self.ions[ia].r[1] * self.reciprocal[2][1] + self.ions[ia].r[2] * self.reciprocal[2][2]
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
	    print "coordinates already in Cartesian"
	    raise ValueError('in Config.dirkar ')
	for ia in range( len ( self.ions ) ) :
	    v1 = self.ions[ia].r[0] * self.direct[0][0] + self.ions[ia].r[1] * self.direct[1][0] + self.ions[ia].r[2] * self.direct[2][0]
	    v2 = self.ions[ia].r[0] * self.direct[0][1] + self.ions[ia].r[1] * self.direct[1][1] + self.ions[ia].r[2] * self.direct[2][1]
	    v3 = self.ions[ia].r[0] * self.direct[0][2] + self.ions[ia].r[1] * self.direct[1][2] + self.ions[ia].r[2] * self.direct[2][2]
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
        except IOError,e:
            print "ERROR : ",e
            print "#"+80*"="
            sys.exit()
        else:
            lines       = f.readlines()

        self.natm   = int(lines[0])
        self.system = lines[1].strip()

       	u           = [ float(item) for item in lines[2].split() ]
        v           = [ float(item) for item in lines[3].split() ]
        w           = [ float(item) for item in lines[4].split() ]
        assert len(u)==3 , "line 3 of POSFF should have only 3 elements"
        assert len(v)==3 , "line 4 of POSFF should have only 3 elements"
        assert len(w)==3 , "line 5 of POSFF should have only 3 elements"

        self.ntype        = int(lines[5])
        self.types        = [ item for item in lines[6].split() ]
        self.natmpertype  = [ item for item in lines[7].split() ]
        self.coord_format = lines[8].strip()

        for i in range(9, 9+self.natm):
            line_list = [ item for item in lines[i].split()]
            self.ions.append (  Ion ( index_ion=i-8, type=line_list[0], pos = [ float(line_list[1]), float(line_list[2]) , float(line_list[3]) ] ) )

        self.volume(u,v,w)
        self.param (u,v,w)
        self.recip (u,v,w)

        f.close()


    # =====================================================================================
    def read_POSCAR(self,filename):

        try:
            f = open(filename)
        except IOError,e:
            print "ERROR : ",e
            print "#"+80*"="
            sys.exit()
        else:
            lines       = f.readlines()

	self.system     = lines[0].strip()
        coeff = float(lines[1])

        u           = [ float(item) for item in lines[2].split() ]
        v           = [ float(item) for item in lines[3].split() ]
        w           = [ float(item) for item in lines[4].split() ]
        assert len(u)==3 , "line 3 of POSFF should have only 3 elements"
        assert len(v)==3 , "line 4 of POSFF should have only 3 elements"
        assert len(w)==3 , "line 5 of POSFF should have only 3 elements"

        self.types        = [ item for item in lines[5].split() ]
#       print self.types
        self.natmpertype  = [ int(item) for item in lines[6].split() ]
        self.ntype        = len(self.natmpertype)
        self.coord_format = lines[7].strip()
	if ( self.coord_format == 'S' or self.coord_format == 'Selective' or self.coord_format == 'Selective dynamics') :
	    self.dyn = True
	    self.coord_format = lines[8].strip()
	    start = 9
	else:
            self.dyn = False
	    start = 8
	self.natm         = sum(self.natmpertype)

        for i in range(start, start+self.natm):
            line_list = [ item for item in lines[i].split()]
            #print line_list
	    if self.dyn :
	        self.ions.append (  Ion ( index_ion=i-8, pos = [ float(line_list[0]), float(line_list[1]) , float(line_list[2]) ] , allowedmove = [ line_list[3] , line_list[4], line_list[5] ] ) )
            else :
	        self.ions.append (  Ion ( index_ion=i-8, pos = [ float(line_list[0]), float(line_list[1]) , float(line_list[2]) ] ) )

	
        self.typeinfo_init()		

        self.lattice_reinit(u,v,w,coeff)
        self.get_volume(u,v,w)
        self.get_lattice_parameters (u,v,w)
        self.get_reciprocal_basis (u,v,w)

	if self.coord_format in self.cartesian_allowed :	
            print "input coordinates in cartesian: apply kardir"
            self.kardir ()

        f.close()

   # =====================================================================================
    def read_multipl_POSCAR(self,f):

        line = f.readline()
        self.system     = line.strip()
        line = f.readline()
        coeff = float(line)

        line = f.readline()
        u           = [ float(item) for item in line.split() ]
        line = f.readline()
        v           = [ float(item) for item in line.split() ]
        line = f.readline()
        w           = [ float(item) for item in line.split() ]
        assert len(u)==3 , "line 3 of POSFF should have only 3 elements"
        assert len(v)==3 , "line 4 of POSFF should have only 3 elements"
        assert len(w)==3 , "line 5 of POSFF should have only 3 elements"

        line = f.readline()
        self.types        = [ item for item in line.split() ]
        line = f.readline()
        self.natmpertype  = [ int(item) for item in line.split() ]
        self.ntype        = len(self.natmpertype)
        line = f.readline()
        self.coord_format = line.strip()
        if ( self.coord_format == 'S' or self.coord_format == 'Selective' or self.coord_format == 'Selective dynamics') :
            self.dyn = True
            line = f.readline()
            self.coord_format = line.strip()
            start = 9
        else:
            self.dyn = False
            start = 8
        self.natm         = sum(self.natmpertype)

        for i in range(start, start+self.natm):
            line = f.readline()
            line_list = [ item for item in line.split()]
            #print line_list
            if self.dyn :
                self.ions.append (  Ion ( index_ion=i-8, pos = [ float(line_list[0]), float(line_list[1]) , float(line_list[2]) ] , allowedmove = [ line_list[3] , line_list[4], line_list[5] ] ) )
            else :
                self.ions.append (  Ion ( index_ion=i-8, pos = [ float(line_list[0]), float(line_list[1]) , float(line_list[2]) ] ) )


        self.typeinfo_init()

        self.lattice_reinit(u,v,w,coeff)

        if self.coord_format in self.cartesian_allowed :
            print "input coordinates in cartesian: apply kardir"
            self.kardir ()


    
    # =====================================================================================
    def write_POSCAR(self,filename):

        f = open(filename,'w+')
	print >> f , self.system
	print >> f , self.coeff
	for i in range(3):	
	    print >> f , float_format*len(self.direct[i]) % tuple(self.direct[i])
	print >> f , '%4s'*self.ntype %  tuple(self.types)
	print >> f , '%4i'*self.ntype %  tuple(self.natmpertype)
	if self.dyn :
            print >> f, 'S'
	    print >> f, self.coord_format
	else:
            print >> f, self.coord_format

	for ia in range(self.natm):
            if self.dyn :
	        form = 3*float_format+3*"%2s"
		print >> f , form % ( self.ions[ia].r[0] , self.ions[ia].r[1], self.ions[ia].r[2] , self.ions[ia].allowedmove[0],  self.ions[ia].allowedmove[1],  self.ions[ia].allowedmove[2] )
            else:
	        form = 3*float_format
        	print >> f , form % ( self.ions[ia].r[0] , self.ions[ia].r[1], self.ions[ia].r[2] )
	f.close()
				
    # =====================================================================================
    def write_POSFF(self,filename):

        f = open(filename,'w+')
        print >> f , self.natm
	print >> f , self.system
	for i in range(3):	
	    print >> f , float_format*len(self.direct[i]) % tuple(self.direct[i])
        print >> f , self.ntype
	print >> f , '%4s'*self.ntype %  tuple(self.types)
	print >> f , '%4i'*self.ntype %  tuple(self.natmpertype)
        print >> f, self.coord_format
	for ia in range(self.natm):
	        form = "%4s"+3*float_format
        	print >> f , form % ( self.ions[ia].type, self.ions[ia].r[0] , self.ions[ia].r[1], self.ions[ia].r[2] )
	f.close()
				

    # =====================================================================================
    def write_TRAJFF(self,filename,option):
    
        f = open(filename,option)
	print >> f , self.natm
        print >> f , self.system
        for i in range(3):
            print >> f , float_format*len(self.direct[i]) % tuple(self.direct[i])
	print >> f , self.ntype
        print >> f , '%4s'*self.ntype %  tuple(self.types)
        print >> f , '%4i'*self.ntype %  tuple(self.natmpertype)
        print >> f, self.coord_format

        for ia in range(self.natm):
	    form = "%4s"+3*float_format
            print >> f , form % ( self.ions[ia].type , self.ions[ia].r[0] , self.ions[ia].r[1], self.ions[ia].r[2] )


# ============================
#        testing unit  :
# ============================
if __name__ == '__main__':

    u = [ 10.0,  0.,  0. ]
    v = [  0.0, 10.,  0. ]
    w = [  0.0,  0., 10. ]
    
    ions=[]
    conf_test=Config(ions=ions,u=u,v=v,w=w)

    print conf_test.print_direct_basis(None)

