#helix axis analysis
#identify helix SSE
#retrieve axis trace for helices
#polynormial fit
#k/t of fitted curve

from math import *
from copy import *
from cacn import *
from sys import *
from rax import *
#from numpy import polyfit
from pnfit import *
from dx import rise
from hlx import *
import numpy as np
import warnings

def axfit(pdb,ax,forder,m,pdbh,ris,dssplocation,hfile,hha,resi):
#  print '#axfit rise',ris
  patha=''

  na=20
  #set ploynormial fit order: 2,3 or 5
#  forder=3

#read helix segments
  ht=helixl(m,pdb,pdbh,dssplocation,hfile,hha,resi)
  

  ds=ris
  nds=len(ds)

#  print '#ht'
#  for a in ht:
#    print a

  #read axis trace position (residue average)
  hxt=raxt(ax,ht)
#  print 'hxt',hxt
  nht=len(ht)
  hfits=[]
  hakts=[]
  r3=[0,1,2]
  warnings.simplefilter('ignore', np.RankWarning)
  for i in range(nht):
#    print ht[i]
    dat=[[],[],[]]
    for a in hxt[i]:
#      print a
      
  #polynormial fit
      for j in r3:
        dat[j].append(a[j])
    nd=len(dat[0])
#    for j in xrange(nd):
#      print dat[0][j],dat[1][j],dat[2][j]
        
    t=range(nd)
  #  print t
    ft=[]
    for j in r3:
      ft1=np.polyfit(t,dat[j],forder)
      ft2=list(ft1)
      ft2.reverse()
#      print ft1
#      print ft2
      ft.append(ft2)
#    print ft
    hfits.append(ft)
  #calculate k,t
    if forder==2:
      hakt=pfit2(ft,nd)
    elif forder==3:
      hakt=pfit3(ft,nd)
    elif forder==5:
      hakt=pfit(ft,nd)
    else:
      print ('use fit order 2,3 or 5')
      break
    hakts.append(hakt)

    
#  print pdb,'done'
  data=[ht,hfits,hakts,hxt]
  return(data)
