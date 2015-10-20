#!/usr/bin/python

import random
from config import Config 
from config import Ion 
from optparse   import OptionParser
import sys

atype=[]
itype=[]
u=[]
v=[]
w=[]

# options parsing 
usage = "usage: %prog --help [-h][-q][-n][-t][-i][-a][-b][-c]"
parser = OptionParser(usage,version="%prog 0.1")
parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,help="don't print status messages to stdout")
# nions
parser.add_option("-n", "--nions", type="int", dest="nions", help="number of ions in the box")
# atype
parser.add_option("-t", "--atype", type='string', dest="satype", help="list of ion types name ( same order that itype)")
# itype
parser.add_option("-i", "--itype", type="string", dest="sitype", help="list of number of ions of each type (same order than atype)")
# u 
parser.add_option("-a", "--acell", type="string", dest="su", help="a vector of the box ")
# v 
parser.add_option("-b", "--bcell", type="string", dest="sv", help="b vector of the box ")
# w
parser.add_option("-c", "--ccell", type="string", dest="sw", help="c vector of the box ")
# cubic
parser.add_option("--cubic", action="store_true", dest="cubic", default=False,help="cubic cell")
# cell
parser.add_option("--cell", type="float", dest="cell", default=None,help="cubic cell parameters")
(options, args ) = parser.parse_args()


if not options.nions :
    parser.error('Number of ions needed --nions')

if not options.satype :
    parser.error('List of ion types name needed --atype')

if not options.sitype :
    parser.error('List of number of ions per type needed --itype')

if options.cubic :
    if not options.cell :
         parser.error('Cubic cell parameters needed --cell')
else:
    if not options.su or not options.sv or not options.sw  :
        parser.error('cell parameters needed --acell --bcell --ccell')


if options.verbose :
    print 60*"="
    print "random_struct in box"
    print "author : fmv"
    print "email  : filipe.manuel.vasconcelos@gmail.com"
    print "email  : filipe.vasconcelos@cea.fr"
    print "date   : mars 2013 "
    print 60*"="

atype=options.satype.split()
itype=[ int(item) for item in options.sitype.split() ]
if options.cubic :
    u    =[ options.cell, 0.0         , 0.0          ]
    v    =[ 0.0         , options.cell, 0.0          ]
    w    =[ 0.0         , 0.0         , options.cell ]
else  :
    u    =[ float(item) for item in options.su.split() ]
    v    =[ float(item) for item in options.sv.split() ]
    w    =[ float(item) for item in options.sw.split() ]

ions=[]	
conf = Config(ions=ions,u=u,v=v,w=w,system='randomPY',natm=options.nions,ntype=len(atype),types=atype, natmpertype = itype, coord_format='Direct')

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
print "POSFF.randomPY generated"



