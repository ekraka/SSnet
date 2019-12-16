#Read DSSP SSE assignment

from math import *
#from readv import *
#from dict import *
from copy import *
#from toolv6h import *

def readssp(proteins):
#Read sequence code
  pdb=proteins
  pathd=''
  protlist=proteins
  y=pathd+proteins+'.dssp'
  res=[[],[],[],[]]#a,b,3-10,pi
  chs=[]

#segment of residues

  f1=open(y,'r')
  fa=f1.readlines()
  f1.close()
  n=len(fa)

  f2=open(pdb+'.ds','w')
  for i in xrange(n):
    if fa[i][:12]=='  #  RESIDUE':
      nst=i+1
      break
  H,E,G,P,B,T,S,U=[],[],[],[],[],[],[],[]
  total=0

  resd=[]
  res=[]
  for i in xrange(nst,n):
#    print 'nst',nst
    a=fa[i]
#    if a[13:15]=='!*':
#      break 
#    lx=fa[i].split()

#    if lx[1]=='!*':
#      break
    b=a[16]
    c=a[13]
    resd.append(b)
    res.append(c)
    f2.write(b+' '+c+'\n')
  f2.close()
#  print res
#  print resd  
  return([res,resd])



  
  

