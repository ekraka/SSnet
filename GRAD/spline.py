#cubic spline fit

from math import *
from utils import *
from copy import *
from sys import *
from skt import *
from figure import *

#spline fit from CA
def cubspline(p,pdc):
#  p=argv[1]
  n=len(pdc)
  max=4.
  eps=1e-15
  eps2=1e-100
  maxt=1e10
  xart=[]
  xart2=[]
  scs2=[]
  dx=.0001
  iter=300
  cbond=3.796484031

  x=[]
  y=[]
  z=[]
  b=[]
  a=[[],[],[]]
  c=[[],[],[]]
  d=[[],[],[]]
  for i in range(n):
    x.append(pdc[i][0])
    y.append(pdc[i][1])
    z.append(pdc[i][2])
  crd=[x,y,z]

  for i in range(3):
    b.append(splcoef(crd[i]))
  #print n
  for i in range(3):
    for j in range(n-2):
  #    print n
      d[i].append(crd[i][j])
      s=crd[i][j+1]-d[i][j]
      a[i].append((b[i][j+1]-b[i][j])/3)
      c[i].append(s-(b[i][j+1]+2.*b[i][j])/3.)
    d[i].append(crd[i][n-2])
    s=crd[i][n-1]-d[i][n-2]
  #  print b[i][n-2]
    a[i].append(-b[i][n-2]/3.)
    c[i].append(s-b[i][n-2]*2/3.)
  dcf=[d,c,b,a]
#  for a in dcf:
#    print 'coef'
#    for b in a:
#      for c in b:
#        print c,
#      print
  sovab=[]
  scs=[]
  skts=sktd2(dcf)
  skts.insert(0,[0.,0.,0.])
  #skts.append([0.,0.,0.])
#  print 'skt'
#  for a in skts:
#    print a
  sa=[0.]
  s=0
  for a in skts:
    s+=a[0]
    sa.append(s)
#  f2=open(p+'.skt','w')
  ns=len(sa)
  i=0
#  for a in skts:
#    i+=1
#    f2.write('%10.4f'%a[0]+'%10.4f'%sa[i]+'%10.4f'%a[1]+'%10.4f'%a[2]+'\n')
#  f2.close()
  return([dcf,skts])

#spline fit from axis traces  
def splinec(p,crdn,seq,resi,patho):
#  p=argv[1]
#  pdb=readcrd(p+'.pdb')
  n=len(crdn)
#  print n
#  print '#spline2 crd',n
#  for a in crdn:
#    print a

  #print dat
  #print len(dat[0])
  max=4.
  eps=1e-15
  eps2=1e-100
  maxt=1e10
  xart=[]
  xart2=[]
  scs2=[]
  dx=.0001
  iter=300
  cbond=3.796484031

  x=[]
  y=[]
  z=[]
  b=[]
  a=[[],[],[]]
  c=[[],[],[]]
  d=[[],[],[]]
  for i in range(n):
    x.append(crdn[i][0])
    y.append(crdn[i][1])
    z.append(crdn[i][2])
  crd=[x,y,z]

  for i in range(3):
    b.append(splcoef(crd[i]))
  #print n
  for i in range(3):
    for j in range(n-2):
  #    print n
      d[i].append(crd[i][j])
      s=crd[i][j+1]-d[i][j]
      a[i].append((b[i][j+1]-b[i][j])/3)
      c[i].append(s-(b[i][j+1]+2.*b[i][j])/3.)
    d[i].append(crd[i][n-2])
    s=crd[i][n-1]-d[i][n-2]
  #  print b[i][n-2]
    a[i].append(-b[i][n-2]/3.)
    c[i].append(s-b[i][n-2]*2/3.)
  dcf=[d,c,b,a]
#  for a in dcf:
#    print 'coef'
#    for b in a:
#      for c in b:
#        print c,
#      print
  sovab=[]
  scs=[]
  skts=sktd2(dcf)
  skts.insert(0,[0.,0.,0.])
  #skts.append([0.,0.,0.])
  
  # f3=open(p+'.axf','w')
  # print 'skt'
  # for a in skts:
  #  print a
  #  for b in a:
  #    f3.write('%10.4f'%b)
  #  f3.write('\n')
  # f3.close()
  
  r3=range(3)
  skts2=[[],[],[]]
  for a in skts:
    for j in r3:
      skts2[j].append(a[j])
  
  sktsc=sktdc(dcf)
  sktsp=sktsc[0]
  crdsp=sktsc[1]
#  print 'crdsp',len(crdsp)
  crdsp.insert(0,crdn[0])
#  print 'ax spline',len(sktsp[0]),len(crdsp),len(skts)
  f2=open(patho+p+'.axsp','w')
  ns=len(sktsp[1])
#  print ns
  [sr,sa,sk,st]=sktsp
#  ns=len(sa)
#  print ns,len(sk),len(crd[0]),len(crdsp)
#  print crdn[0]
  i=0
#  for a in skts:
#  crdsp.insert(0,[x[0],y[0],z[0]])
#  
  
  
  f2.write('#     si         s         k         t        x         y         z        kCA       tCA    Res  #\n')
  for i in range(ns):
#    i+=1
    crx=crdsp[i]
#    print crx
    f2.write('%10.4f'%sr[i]+'%10.4f'%sa[i]+'%10.4f'%sk[i]+'%10.4f'%st[i]\
      +'%10.4f'%crx[0]+'%10.4f'%crx[1]+'%10.4f'%crx[2])
    if i%20==0:
      j=int(i/20)
      f2.write('%10.4f'%sk[i]+'%10.4f'%st[i]+'%4s'%seq[j+1]+'%4s'%resi[j+1])
#      f3.write('%10.4f'%sr[i]+'%10.4f'%sa[i]+'%10.4f'%sk[i]+'%10.4f'%st[i]+'\n')
    f2.write('\n')
  f2.close()
#  f3.close()
#  print sa[-1]
#  wsgp(p,sa[-1])
#  print len(crdsp)

  return([dcf,sktsp,crdsp,skts2])


