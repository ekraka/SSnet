#read 2' structure elements definetion
from sys import *
def ssedf(pdb,ch,pathp):
  name=pdb+ch
  f1=open(pathp+pdb+'.pdb')
  f2=open(name+'p.pdb')
  fa=f1.readlines()
  f1.close()
  fb=f2.readlines()
  f2.close()
  na=len(fa)
  nb=len(fb)
  ch=fb[1][21]
  sh,ss=[],[]
  H,S=[],[]
  for a in fa:
    la=a.split()
    if la[0]=='HELIX':
  #    ch1=a[19:20]
      if a[19:20]==ch:
#        H.append([la[5],la[8]])
        H.append([a[21:26],a[33:38]])
        sh.append(a)
  #    print la[5],la[8]
    if la[0]=='SHEET':
  #    print la
      if a[21:22]==ch:
        S.append([a[22:27],a[33:38]])
        ss.append(a)
  #      print la[6],la[9]
    if la[0]=='ATOM':
      break

  pd=[]
  CA=[]
  resi=[]
  for a in fb:
    if len(a)<2:
      fb.remove(a)
  #c=0
  #end=0
#  check=fb[0][22:26]
  n=[]
  m=[]
#  n.append(int(check))
  for a in fb:
#    la=a.split()
#    if la[0]=='ATOM':
    if a[13:15]=='CA':
      pd.append(a)
      resi.append(a[22:27])


#  print H
  nh=len(H)
    
  H2=[]
#  sh2=[]
  for a in H:
#    print a
#change id to res # in pdb by +1
    H2.append([resi.index(a[0])+1,resi.index(a[1])+1])
#  print 'H2',H2 
  S2=[]
  ss2=[]

  for a in S:
    S2.append([resi.index(a[0])+1,resi.index(a[1])+1])
#  print 'S2',S2
  ns=len(S)

  AA={'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E',\
      'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',\
      'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W',\
      'TYR':'Y', 'VAL':'V'}
  
#read crd of CA
  crd=[]
  seq=[]
  resi=[]
  for a in pd:
    crd.append([float(a[30:38]),float(a[38:46]),float(a[46:54])])
    res=AA.get(a[17:20])
    resi.append(a[22:26])
    if res!=None:
      seq.append(res)
    else:
      seq.append('X')
    
  #print sse
  #print pd
  f2=open(name+'s.pdb','w')
  for i in range(nh):
#    print H2[i][0],H2[i][1]
    b=sh[i][:20]+'%5d'%H2[i][0]+' '+sh[i][26:32]+'%5d'%H2[i][1]+''+sh[i][38:]
    f2.write(b)
  for i in range(ns):
    b=ss[i][:22]+'%4d'%S2[i][0]+ss[i][26:33]+'%4d'%S2[i][1]+ss[i][37:]
    f2.write(b)
  for a in pd:
  #  print a
    f2.write(a)
#  print pdb, ' new pdb created'
  f2.close()
  return([H2,crd,seq,resi])
