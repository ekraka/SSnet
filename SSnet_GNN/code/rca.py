#locate CA coordinates
from math import *
def readca(pdb,path):
#  print 'read data'

  f1=open(path+pdb+'.pdb')
  fa=f1.readlines()
  ln=len(fa)
  pd=[]
  N=[]
  CA=[]
  C=[]
  c=0
  end=0
  for a in fa:
    la=a.split()
    if la[0]=='ATOM':
      st=c
      break
    c+=1
  c=0
  for a in fa:
    la=a.split()
#      ln=int(la[5].strip())
    if la[0]=='TER':
      end=c
      break
    c+=1
  if end==0:
    end=ln
#  print st,end
#  print 'c', c
  for i in range(st,end):
    a=fa[i]
    la=a.split()
    if la[2]=='CA':
#      print a[30:38],a[38:46],a[46:54]
      CA.append([float(a[30:38]),float(a[38:46]),float(a[46:54])])

  nl=len(CA)
  return(CA)

