#get helix segments

import os
import os.path
from dscode import readssp

def helixl(m,pdb,pdbh,dssplocation,hfile,hha,resi):
  if m==1:
#  hlx=hpdb(pdb)
    hlx=pdbh
  elif m==2:
    hlx=hdssp(pdb,dssplocation)
  elif m==3:
    hlx=hcustm(hfile,resi)
#     if hlx=='':
#       print 'no ',hfile,' found. Use defintion in pdb'
#       hlx=pdbh
  elif m==0:
#    hlx=hapsa(pdb)
    hlx=hha
  else:
#    hlx=[]
    print ('Note! no valid helix definition given, use defintion in pdb')
    hlx=pdbh
#  print 'hlx',hlx
  return(hlx)
 
#read SSE definition from pdb
def hpdb(pdb):
  paths=''
  f1=open(paths+pdbs+'.pdb','r')
  fb=f1.readlines()
  f1.close()
  nfb=len(fb)
  H,B=[],[]
  for i in xrange(nfb):
    a1=fb[i].split()
    if a1[0]=='ATOM':
      break
    if a1[0]=='HELIX':
      H.append([a1[5],a1[8]])
    elif a1[0]=='SHEET':
      B.append([a1[6],a1[9]])
  #print H, B
#  print 'H'
#  for a in H:
#    print a
#  print 'B'
#  for a in B:
#    print a
  ht=[]
  if H!=[]:
    for a in H:
      ht.append([a[0]+1,a[-1]+1])
  return(ht)
  
#read dssp helix
def hdssp(pdb,pathd):
#read sse defintion from dssp, store in he
#  cmd='../tool/dsspe'+' '+name+'p.pdb '+name+'.dssp'
  cmd=pathd+' '+pdb+'p.pdb '+pdb+'.dssp'
#  else:
#  cmd=path+' -ch '+chain+' '+pdb+'p.pdb '+'/dev/stdout'
#  print cmd
  p=os.popen(cmd)
#  f1=open(pathd+pdb+'.ds')
#  fd=f1.readlines()
#  f1.close()
  dsspc=readssp(pdb)
#  print 'dssp read',dsspc
#  nfd=len(fd)
#  D=[]
#  for a in fd:
#    D.append(a[0])
  #if nfd!=na3:
  #  print '#res diff A&D',pdb,na3,nfd
  D=dsspc[1]
#  print D
  dssp=['H','E','G','I','B','S','T',' ']
  dsdc={'H':0,'E':1,'G':2,'I':3,'B':4,'S':5,'T':6,' ':7}
  di=[[],[],[],[],[],[],[],[]]
  ks=[[],[],[],[],[],[],[],[]]
  ks0=[[],[],[],[]]
  kst=[[],[],[],[]]
  nfd=len(D)
  for i in xrange(nfd):
    di[dsdc[D[i]]].append(i)
  #print di
  ksm=[]
  ndi=len(dssp)

  hxd=di[0]
#  print hxd
  nda=len(hxd)
  he=[]
  if nda>0:
    d0=hxd[0]
    bd=[d0]
    for j in xrange(1,nda):
      if hxd[j]-d0==1:
        bd.append(hxd[j])
        d0=hxd[j]
      else:
        d0=hxd[j]
        he.append(bd)
        bd=[d0]
    if bd!=[]:
      he.append(bd)
#   print 'helix dssp'
#   for a in he:
#    print a
  ht=[]
  if he!=[]:
    for a in he:
      ht.append([a[0]+1,a[-1]+1])
  return(ht)

#read customized helix definition
def hcustm(f,resi):
  pathf=''
  fexist=os.path.isfile(pathf+f)
  if fexist==False:
    print ('no ',f,' found')
    return([])
  f1=open(pathf+f,'r')
  fa=f1.readlines()
  f1.close()
  hx=[]
  hx2=[]
  for a in fa:
    b=a.split()
    if len(b)==2:
      c=[]
      c3=[]
      for c1 in b:
        c.append(int(c1))
        c1='%4s'%c1
        c2=resi.index(c1)+1
        c3.append(c2)
      hx.append(c)
      hx2.append(c3)
#   print hx
#   print hx2
  return(hx2)
  
def dssprer(pdb,pathd):
#read sse defintion from dssp, store in he
#  cmd='../tool/dsspe'+' '+name+'p.pdb '+name+'.dssp'
  cmd=pathd+' '+pdb+'p.pdb '+pdb+'.dssp'
#  else:
#  cmd=path+' -ch '+chain+' '+pdb+'p.pdb '+'/dev/stdout'
#  print cmd
  p=os.popen(cmd)
#  f1=open(pathd+pdb+'.ds')
#  fd=f1.readlines()
#  f1.close()
