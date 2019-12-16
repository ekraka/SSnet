#calculate rise per residue

from math import *
from copy import *
#from cacn import *
from sys import *
from tool import *

def rise(ax):
  crd=[]
  crd=ax
#  for a in crd:
#    print a
  nc=len(crd)

  dst=[]
  for i in xrange(nc-1):
    a=crd[i]
    b=crd[i+1]
    dst.append(dist(a,b))

#  print 'distance'
#  f2=open(pdb+'.dis','w')
#  for a in dst:
#    print a
#    f2.write('%10.6s'%a+'\n')
#  f2.close()
  return(dst)

