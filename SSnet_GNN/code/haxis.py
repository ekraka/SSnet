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
from rca import *
import os
from pdb import *
from ssdh import *
from io import StringIO
#import cStringIO
from vech import *
from spline import *
from dl import *
from ha import *
from dx import *
from output import *
from hct import *
from kink import *
from kinkd import *
from hhas import *
from dia import *

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
#    print seq
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
    '''
    #exit()
#    for a in skts:
#      print a
    coefs=splined[0]
    coefs.reverse()
    coef=coefs[:3]
#    print '***coef***'
#    for a in coefs:
#      for b in a:
#        print b
#      print
  #  crd2=splined[2]
  #  print '#crd2'
  #  for a in crd2:
  #    print a

    vtnb=vecs(coef)
#    print '***TNB***'
#    for a in vtnb:
#      for b in a:
#        print b
#      print

    axs=axis(vtnb,crd)
    ax=axs[0]
#    print '***ax***'
#    print ax

    ris=rise(ax)
#    print '***ris***',ris
#    print 'spline fit of axis'
    splineax=splinec(name,ax,seq,resi,patho)
    skta=splineax[1]
    sktc=splineax[3]
  #  m=1
  #  dssplocation='../dssp'
  #  forder=3
#    print 'name',pdb,name
    hfile=name+'.hx'
    diamt=diam(crd,ax,name)
    hha=hhs(pdbh,skta,ris,name,resi,sktc,diamt,axs)
    md=[mi]
    if dual==1:
      md.append(0)
    o2=0
#    print md
    for m in md:
      axft=axfit(name,ax,2,m,pdbh,ris,dssplocation,hfile,hha,resi)
      hd=axft[0]
  #    print axft[0]
      if axft[0]==[]:
        print 'no helix. Skip'
        break
      axft3=axfit(name,ax,3,m,pdbh,ris,dssplocation,hfile,hha,resi)
      axft5=axfit(name,ax,5,m,pdbh,ris,dssplocation,hfile,hha,resi)
      hlxc=gethclass(axft)
  #    print 'hclass',hlxc
      kr=rish(ris,hd,nca)
      kq=kinkq(hd,skta[2],ris,name,sktc)
  #    print kr
      o2+=1
      outputh(name,axft,ris,ax,vtnb,skts,pre0,pdbcrd,hlxc,kr,kq,axft3,axft5,hha,resi,skta[1],dual,o2,patho)
    clean(name)
    '''
    print (name, 'done')
    
if __name__ == "__main__":
  haxis()

