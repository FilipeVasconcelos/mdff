import sys
import numpy as np
from constants import float_format

class Lattice(object):
    """
    Lattice class: simple definitions.
    ------
        lattice related parameters from the direct basis vectors
    
    """

    # intialize the direct basis vector
    def __init__( self, u = [ 1., 0., 0. ], v = [ 0., 1., 0. ], w = [ 0., 0., 1. ], coeff = 1. ):
        """
        Inititialization of the direct basis vectors.
        use method cell_change() to update the cell parameters for
        a new set of direct basis vectors
        """  
	direct                       = [ u, v, w ]
        self.coeff                   = coeff
	self.direct                  = [ [ i*coeff for i in j ] for j in direct ]
        self.reciprocal_updated      = False
        self.lattice_parameters_updated = False 
        self.volume_updated          = False
        self.get_lattice_parameters( u, v, w )
        self.get_volume( u, v, w )
        self.get_reciprocal_basis( u, v, w )

    def lattice_change( self, du, dv, dw ):

        for i in range ( 3 ) :
            self.direct[0][i] += du[i]
        for i in range ( 3 ) :
            self.direct[1][i] += dv[i]
        for i in range ( 3 ) :
            self.direct[2][i] += dw[i]
        self.reciprocal_updated      = False
        self.lattice_parameters_updated = False
        self.volume_updated          = False

    def lattice_reinit( self, u, v, w , coeff ):

        for i in range ( 3 ) :
            self.direct[0][i] = u[i]*coeff
        for i in range ( 3 ) :
            self.direct[1][i] = v[i]*coeff
        for i in range ( 3 ) :
            self.direct[2][i] = w[i]*coeff
        self.reciprocal_updated         = False
        self.lattice_parameters_updated = False
        self.volume_updated             = False
        self.get_lattice_parameters(  self.direct[0],  self.direct[1],  self.direct[2] )
        self.get_volume( self.direct[0],  self.direct[1],  self.direct[2]  )
        self.get_reciprocal_basis( self.direct[0],  self.direct[1],  self.direct[2]  )

    def print_direct_basis( self, filename ):
        for i in range( 3 ):	
            print >> filename , float_format*len(self.direct[i]) % tuple(self.direct[i])

    def print_reciprocal_basis( self, filename ):
        if self.reciprocal_updated :
            for i in range ( 3 ):	
                print >> filename , float_format*len(self.reciprocal[i]) % tuple(self.reciprocal[i])
        else:
            u = self.direct[0]
            v = self.direct[1]
            w = self.direct[2]
            self.get_reciprocal_basis( u, v, w )
            for i in range ( 3 ):
                print >> filename , float_format*len(self.reciprocal[i]) % tuple(self.reciprocal[i])

    def print_lattice_parameters( self, filename ):
        if self.lattice_parameters_updated :   
            print >> filename , ("%4s"+float_format)*3 % ("a      = ",self.lattice_lengths[0]    ," b    = ",self.lattice_lengths[1]   ," c     = ",self.lattice_lengths[2] )
            print >> filename , ("%4s"+float_format)*3 % ("alpha  = ",self.lattice_angles[0]     ," beta = ",self.lattice_angles[1]    ," gamma = ",self.alttice_angles[2] )
        else:
            u = self.direct[0]
            v = self.direct[1]
            w = self.direct[2]
            self.get_lattice_parameters( u, v, w )
            print >> filename , ("%4s"+float_format)*3 % ("a      = ",self.lattice_lengths[0]    ," b    = ",self.lattice_lengths[1]   ," c     = ",self.lattice_lengths[2] )
            print >> filename , ("%4s"+float_format)*3 % ("alpha  = ",self.angles[0]             ," beta = ",self.lattice_angles[1]    ," gamma = ",self.lattice_angles[2] )

    def print_volume(self, filename):
        if self.volume_updated :
            print >> filename, ("%s"+float_format) % ("volume = ",self.volume)
        else:
            u = self.direct[0]
            v = self.direct[1]
            w = self.direct[2]
            self.get_volume( u, v, w )
            print >> filename, ("%s"+float_format) % ("volume = ",self.volume)
            
	
    # get cell parametes a,b,c,alpha,beta,gamma 
    def get_lattice_parameters( self, u, v, w ):

        a = np.linalg.norm( u )
        b = np.linalg.norm( v )
        c = np.linalg.norm( w )

        alpha = w[0] * v[0] + w[1] * v[1] + w[2] * v[2]
        beta  = u[0] * w[0] + u[1] * w[1] + u[2] * w[2]
        gamma = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
        alpha = alpha / ( c * b )
        beta  = beta  / ( a * c )
        gamma = gamma / ( a * b )
        alpha = np.degrees( np.arccos ( alpha ) )
        beta  = np.degrees( np.arccos ( beta  ) )
        gamma = np.degrees( np.arccos ( gamma ) )

        self.lattice_lengths = [ a, b, c ]
        self.lattice_angles = [ alpha, beta, gamma ]
        self.lattice_parameters_updated = True

    def get_reciprocal_basis( self, u, v, w ):

        self.get_volume( u, v, w )
        urec = np.cross( v, w ) / self.volume
        vrec = np.cross( w, u ) / self.volume
        wrec = np.cross( u, v ) / self.volume
        self.reciprocal = [ urec , vrec, wrec ]

        self.reciprocal_basis_updated = True

    def get_volume( self, u, v, w ):
        urec = np.cross( v, w )
        vrec = np.cross( w, u )
        wrec = np.cross( u, v )
        self.volume = urec[0] * u[0] + vrec[0] * v[0] + wrec[0] * w[0]

        wa = self.volume / np.sqrt ( urec[0] * urec[0] + urec[1] * urec[1] + urec[2] * urec[2] )
        wb = self.volume / np.sqrt ( vrec[0] * vrec[0] + vrec[1] * vrec[1] + vrec[2] * vrec[2] )
        wc = self.volume / np.sqrt ( wrec[0] * wrec[0] + wrec[1] * wrec[1] + wrec[2] * wrec[2] )
        self.dist_opposite_faces = [ wa, wb, wc ]

        self.volume_updated = True

# ============================
#        testing unit  :
# ============================
if __name__ == '__main__':

    u = [ 10.0,  0.,  0. ]
    v = [  0.0, 10.,  0. ]
    w = [  0.0,  0., 10. ]
    box = Lattice( u, v, w )

    print __name__
    print 30*'=' 
    print "testing unit of class: Lattice"
    print  box.__doc__
    print 30*'=' 
    print 30*'-' 
    print "First test cubic cell: " 
    print 30*'-' 
    print "input direct basis :"
    box.print_direct_basis ( filename=None ) # standard output filename = None
    print
    print "cell parameters : "
    box.print_lattice_parameters ( filename=None ) 
    print
    box.print_volume( filename=None )
    print
    print "reciprocal basis :"
    box.print_reciprocal_basis ( filename=None ) # standard output filename = None
    
    du = [  1.0,  0.05,  0. ]
    dv = [  0.0,  1.,  0.01 ]
    dw = [  0.0,  0.01 ,  1. ]
    box.cell_change( du, dv, dw )
    print "input direct basis :"
    box.print_direct_basis ( filename=None ) # standard output filename = None
    print
    print "cell parameters : "
    box.print_lattice_parameters ( filename=None ) 
    print
    box.print_volume( filename=None )
    print
    print "reciprocal basis :"
    box.print_reciprocal_basis ( filename=None ) # standard output filename = None
    
    
