import sys
import numpy as np
from constants import float1, float3, Sep,sep

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
        use method lattice_change() to update the cell parameters for
        a new set of direct basis vectors
        """  
        direct                          = [ u, v, w ]
        self.coeff                      = coeff
        self.direct                     = [ [ i*coeff for i in j ] for j in direct ]
        self.lattice_parameters( self.direct )
        self.volumeANDoppositefaces( self.direct )
        self.reciprocal_basis( self.direct )

    def lattice_change(self, du, dv, dw ):

        for i in range (3) :
            self.direct[0][i] += du[i]
            self.direct[1][i] += dv[i]
            self.direct[2][i] += dw[i]
        self.lattice_parameters( self.direct )
        self.volumeANDoppositefaces( self.direct )
        self.reciprocal_basis( self.direct )

    def lattice_reinit( self, u, v, w , coeff ):
        for i in range (3) :
            self.direct[0][i] = u[i]*coeff
            self.direct[1][i] = v[i]*coeff
            self.direct[2][i] = w[i]*coeff
        self.lattice_parameters( self.direct )
        self.volumeANDoppositefaces( self.direct )
        self.reciprocal_basis( self.direct )

    def print_direct_basis( self, f ):
        if f :
            for i in range(3):
                for j in range(3):
                    f.write(float1.format(self.direct[i][j]))
                f.write('\n')
        else:
            for i in range (3):	
                print(float3.format(self.direct[i]))

    def print_reciprocal_basis( self, f):
        if f :
            for i in range(3):
                for j in range(3):
                    f.write(float1.format(self.reciprocal[i][j]))
                f.write('\n')
            f.write('\n')
        else:
            for i in range(3):	
                print(float3.format(self.reciprocal[i]))

    def print_lattice_parameters( self, filename):
        strout1=["a     =","b     =","c     ="]
        strout2=["alpha =","beta  =","gamma ="]
        for i in range(3):
            if f :
                f.write(strout1[i]+float1.format(self.lattice_lengths[i]))
                f.write(strout2[i]+float1.format(self.lattice_angles[i]))
            else:
                print(strout1[i]+float1.format(self.lattice_lengths[i]))
                print(strout2[i]+float1.format(self.lattice_angles[i]))

    def print_volume(self, f):
        if f :
            f.write("volume = "+float1.format(self.volume))
        else:
            print("volume = "+float1.format(self.volume))
            
	
    # get cell parametes a,b,c,alpha,beta,gamma 
    def lattice_parameters( self, basis ):

        u = basis[0]
        v = basis[1]
        w = basis[2]
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

    def reciprocal_basis( self, basis ):
        u=basis[0];v=basis[1];w=basis[2]
        self.volumeANDoppositefaces( basis )
        urec = np.cross( v, w ) / self.volume
        vrec = np.cross( w, u ) / self.volume
        wrec = np.cross( u, v ) / self.volume
        self.reciprocal = [ urec , vrec, wrec ]

    def volumeANDoppositefaces( self, basis):
        u=basis[0];v=basis[1];w=basis[2]
        urec = np.cross( v, w )
        vrec = np.cross( w, u )
        wrec = np.cross( u, v )
        self.volume = urec[0] * u[0] + vrec[0] * v[0] + wrec[0] * w[0]
        wa = self.volume / np.sqrt ( urec[0] * urec[0] + urec[1] * urec[1] + urec[2] * urec[2] )
        wb = self.volume / np.sqrt ( vrec[0] * vrec[0] + vrec[1] * vrec[1] + vrec[2] * vrec[2] )
        wc = self.volume / np.sqrt ( wrec[0] * wrec[0] + wrec[1] * wrec[1] + wrec[2] * wrec[2] )
        self.dist_opposite_faces = [ wa, wb, wc ]

# ============================
#        testing unit  :
# ============================
if __name__ == '__main__':

    u = [ 10.0,  0.0,  0.0 ]
    v = [  0.0, 10.0,  0.0 ]
    w = [  0.0,  0.0, 10.0 ]
    box = Lattice( u, v, w )

    print(sep) 
    print("testing unit of class: Lattice")
    print(sep) 
    print( box.__doc__)
    print(sep) 
    print("First test cubic cell: " )
    print(sep) 
    print("input direct basis :")
    box.print_direct_basis(None) # standard output filename = None
    print()
    print("lattice parameters : ")
    box.print_lattice_parameters (filename=None ) 
    print()
    box.print_volume( filename=None )
    print()
    print("reciprocal basis :")
    box.print_reciprocal_basis (None) # standard output filename = None
    
    du = [  1.0,  0.05,  0. ]
    dv = [  0.0,  1.,  0.01 ]
    dw = [  0.0,  0.01 ,  1. ]
    box.lattice_change( du, dv, dw )
    print("input direct basis :")
    box.print_direct_basis (None ) # standard output filename = None
    print()
    print("cell parameters : ")
    box.print_lattice_parameters ( filename=None ) 
    print()
    box.print_volume( filename=None )
    print()
    print("reciprocal basis :")
    box.print_reciprocal_basis (None) # standard output filename = None
    
    
