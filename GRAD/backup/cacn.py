from math import *
def readcacn(pdb):
#  print 'read data'

  f1=open(pdb+'.pdb')
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
  print (st,end)
  print ('c', c)
  for i in range(st,end):
    a=fa[i]
    la=a.split()
    if la[2]=='N':
      N.append([float(a[31:38]),float(a[39:46]),float(a[47:54])])
    if la[2]=='CA':
      CA.append([float(a[31:38]),float(a[39:46]),float(a[47:54])])
    if la[2]=='C':
      C.append([float(a[31:38]),float(a[39:46]),float(a[47:54])])
  print ('N,CA,C')
  nl=len(CA)
  for i in range(nl):
    for b in [N,CA,C]:
      print ('%8.3f'%b[i][0],'%8.3f'%b[i][1],'%8.3f'%b[i][2])
  return([N,CA,C])
#    print '[', la[23:27],']'
#    print '[', la[23:27].strip(),']'
#    ln=int(la[23:27].strip())

#  f2=open(pdb+'cacn','w')
#print pls
#  for i in range(len(pls)):
#    if pls[i]!=pls[i-1]:
#      f2.write(pls[i]+'\n')
#    i+=1
#  f2.close()
