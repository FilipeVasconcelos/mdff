import numpy as np

Sep=30*"="
sep=30*'-'
float1       = '{0:16.12f} '
float3       = '{0[0]:>16.12f} {0[1]:>16.12f} {0[2]:>16.12f}'
posff        = '{0:2s} {1:16.12f} {2:16.12f} {3:16.12f}\n'
poscar       = '{0:16.12f} {1:16.12f} {2:16.12f} {3:2s}\n'
config       = '{0:16.12f} {1:16.12f} {2:16.12f}\n'
float_format_config = '%20.9f'

piroot       = np.sqrt(np.pi) 
tpi          = 2. * np.pi
fpi          = 2. * tpi
imag         = 1j
rytoev       = 13.605826                        
hart         = rytoev*2.0                      
bohr         = 0.529177249                    
e_2          = hart * bohr      
