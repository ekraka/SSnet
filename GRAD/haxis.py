############################################
#HAXIS main program
#version 0.0.1a
#need input file:
#line1: #full path where pdb files are located, Tab, option of helix defition (0-2), Tab, path to dssp executable if previous option is 1
#line2...:pdb name, Tab, chain id, Tab, model id
#by Daniel Guo at CATCO at SMU
#Mar. 2013
############################################

import sys
from math import *
from tool import *
import os
from pdb import *
from ssdh import *
from io import StringIO
from spline import *

def haxis(inp):
  pathapsa=''
  pathpdb=''
  patho=''

  #read protein list
  #if len(sys.argv)==1:
  #  print 'protein list file needed:'
  #  print 'python apsa_path/apsa2.py protein_list'
  #  exit()
  finput=inp #sys.argv[1]
  pre=readin(finput)
  #locate given chain/model
#  print pre
  pre0=pre[0]
  pdbs=pre[1]
  pathp,mi=pre0[0],pre0[1][0]
  patho=pre0[2]
  dual=pre0[1][1]
#  if pre0[2]!='':
  dssplocation=''
  if mi==2:
    dssplocation=pre0[3]
    
#  if m==1 and pre0[2]=='':
#    print 'dssp executable path required'
#    exit()
  #print pdbs
  err=[]
  #chs=[]
  for p in pdbs:
    pdb=p[0]
    if len(p)>1:
      if p[1]!=',':
        ch=p[1]
      else:
        ch=''
    else:
      ch=''
    if len(p)==3:
      mo=p[2]
    else:
      mo=''
  #  print 'ch',ch,'mo',mo
    fh = StringIO()
    if mo=='':
      chain=pdbp2(pdb,ch,pathp,fh)
      if chain=='':
        print ('no chain ',ch,' found in ',pdb)
        continue
  #    print 'ch',chain
      if chain==' ':
        chain=''
      pdbcrd=ssedf(pdb,chain,pathp)
      name=pdb+chain
    else:
      print ('NMR model',mo)
      model=pdbm(pdb,ch,mo,pathp,fh)
      if model=='':
        print ('no model ',mo,' found in ',pdb)
        continue
      pdbcrd=ssedf(pdb,model,pathp)
      name=pdb+model
  #  chs.append(chain)
    pdbh=pdbcrd[0]
    crd=pdbcrd[1]
    seq=pdbcrd[2]
    resi=pdbcrd[3]
#    print resi
#    print (seq)
#    print '#crd'
#    for a in crd:
#      print a
    fh.close()
    if mi==2:
      dssprer(name,dssplocation)
#      hdssp(name,dssplocation)
      
#    print 'pdbh'
#    print pdbh
#    print name
    nca=len(crd)
    splined=cubspline(name,crd)
    skts=splined[1]
    f = open(name+'.htxt', 'w')
    for abc in skts:
      f.write(' '.join(list(map(str, abc)))+'\n')
    f.close()


    print (name, 'done')
    
if __name__ == "__main__":
  haxis()

