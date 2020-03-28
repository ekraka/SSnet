#Toolkit of functions
from math import *
import os
import os.path

#contact matrix
def contactm1(crd):
#  for a in crd:
#    print a
  cas=len(crd)
#  print cas
  dis=[]
  i=0
  for i in range(cas):
#    print i
    dr=[]
    for k in range(cas):
#      print i,k
      x=0
      for j in range(3):
        t=crd[k][j]-crd[i][j]
        x+=t*t
      xs=sqrt(x)
#contact identified if xs less than cutoff 10
      if xs<10:   
        dr.append(xs)
      else:
        dr.append(0)
    dis.append(dr)
#  i=0
#  for a in dis:
#    for a1 in a:
#      print str('%4.1f'%a1),
#    print
#    i+=1              
  return(dis)

#read chain and model input
def pname(argv):
  m=''
  ch=''
  argt=''
  np=(len(argv)-2)/2
  if np >= 1:
    for i in range(np):
      if argv[2+i*2]=='-m':
        m=argv[3+i*2]
        argt+=' -mo '+m+' '
      elif argv[2+i*2]=='-c':
        ch=argv[3+i*2]
        argt+=' -ch '+ch+' '
  mc=m+ch
  return([argt,mc,m,ch])

#identify program and pdb file location
def patha(argv):
  pathapsa=''
  pathpdb=''
  path=argv[0]
  for i in range(len(path)-1,-1,-1):
    if path[i]=='/':
      pathapsa=path[:i+1]+'apsa'
      break

  pdbf=argv[1]
  for i in range(len(pdbf)-1,-1,-1):
    if pdbf[i]=='/':
      pathpdb=pdbf[:i+1]
      pdb=pdbf[i+1:]
      break
  if pathpdb=='':
    pdb=pdbf
  return(pathapsa,pathpdb,pdb)

#process protein list
def plist(plsts):
  f1=open(plsts,'r')
  fa=f1.readlines()
  pdbs=[]
  pathpdb=''
  path=fa[0].split()
  if len(path)==1:
    pathapsa=path[0]+'apsa'
  elif len(path)==2:
    pathapsa=path[0]+'apsa'
    pathpdb=path[1]
  for line in fa[1:]:
    if line!='\n':
      la=line.split()
      pdbs.append(la)
  f1.close()
  return(pathapsa,pathpdb,pdbs)

#identify chain, model from pdb list file
def papre(p):
  ch,mo='',''
  pdb=p[0]
  np=len(p)
  argt=''
  if np>=2 and p[1]!=',':
    ch=p[1]
    argt+=' -ch '+ch+' '
  if np==3:
    mo=p[2]
    argt+=' -mo '+mo+' '
  mc=mo+ch
  name=pdb+mc
  return(pdb,argt,name)

#coordinate transfer
def ctran(C,p):
  new=[]
  for i in range(3):
    w=0
    for j in range(3):
      w+=C[i][j]*p[j]
    new.append(w)
  return(new)

#caculate euler angles in zxz formulation
def euzxzs(Z):
#  eu=[]
#  print Z[1]
  if abs(Z[1]-1)>1.e-15:
#    print Z[1],'not 1'
    dv=sqrt(1-Z[1]**2)
    alpha=acos(-Z[0]/dv)
    beta=acos(Z[1])
    gama=acos(Z[2]/dv)
    if Z[3]<0:
      alpha=-alpha
    if Z[4]<0:
      gama=-gama 
  else:
#    dv=sqrt(1-Z[1]**2)
#    print '=0'
    alpha=0
#    beta=acos(Z[1])
    beta=0
    gama=0
  eu=[alpha,beta,gama]
#  eun=[-alpha,-beta,-gama]
#  eun=[2*pi-alpha,2*pi-beta,2*pi-gama]
#  eu=[eup,eun]
  return(eu)

#determine two angles (in plane and out of plane angle)
def plane2ang(V):
  beta=asin(V[2])
#  if V[0]>=0:
#    beta=be
#  else:
#    if be>=0:
#      beta=pi-be
#    else:
#      beta=-pi-be
  if abs(V[0]-0)>e-15:
    alpha=atan2(V[1],V[0])
  else:
    if V[1]>=0:
      alpha=pi/2
    else:
      alpha=-pi/2
  return([alpha,beta])

#read spline coefficients
def readcoef(fname):
  f1=open(fname,'r')
  i=0
  j=0
  fa=f1.readlines()
  n=len(fa)
  for i in range(n):
    if fa[i]=='# Spline fit curve equation coefficients: XYZ\n':
      n1=i+2
  #    print n1
    if fa[i][0:2]=='#Y':
      n0=i-n1
    if fa[i]=='## double precision a,k,t\n':
      n2=i+2
  abc=[]
  cfx,afx,bfx=[],[],[]
#  print 'n0',n0
  for i in range(n0):
    la=fa[n1+i].split()
    cfx.append(float(la[0]))
    afx.append(float(la[2]))
    bfx.append(float(la[1]))
  cfy,afy,bfy=[],[],[]
  for i in range(n0):
    la=fa[n1+n0+1+i].split()
    cfy.append(float(la[0]))
    afy.append(float(la[2]))
    bfy.append(float(la[1]))
  cfz,afz,bfz=[],[],[]
  for i in range(n0):
    la=fa[n1+2*n0+2+i].split()
    cfz.append(float(la[0]))
    afz.append(float(la[2]))
    bfz.append(float(la[1]))
    
#  print 'spline coefficient data read'
  f1.close()
  return([[afx,afy,afz],[bfx,bfy,bfz],[cfx,cfy,cfz]])

#read spline coefficients
def readcoefv(t):
#  f1=open(fname,'r')
  i=0
  j=0
#  fa=f1.readlines()
  fa=t
  n=len(fa)
  for i in range(n):
    if fa[i]=='# Spline fit curve equation coefficients: XYZ\n':
      n1=i+2
  #    print n1
    if fa[i][0:2]=='#Y':
      n0=i-n1
    if fa[i]=='## double precision a,k,t\n':
      n2=i+2
  abc=[]
  cfx,afx,bfx=[],[],[]
#  print 'n0',n0
  for i in range(n0):
    la=fa[n1+i].split()
    cfx.append(float(la[0]))
    afx.append(float(la[2]))
    bfx.append(float(la[1]))
  cfy,afy,bfy=[],[],[]
  for i in range(n0):
    la=fa[n1+n0+1+i].split()
    cfy.append(float(la[0]))
    afy.append(float(la[2]))
    bfy.append(float(la[1]))
  cfz,afz,bfz=[],[],[]
  for i in range(n0):
    la=fa[n1+2*n0+2+i].split()
    cfz.append(float(la[0]))
    afz.append(float(la[2]))
    bfz.append(float(la[1]))
    
#  print 'spline coefficient data read'
#  f1.close()
  return([[afx,afy,afz],[bfx,bfy,bfz],[cfx,cfy,cfz]])

#vector cross product
def multcros(ma,mb):
  m1=ma[1]*mb[2]-mb[1]*ma[2]
  m2=ma[2]*mb[0]-ma[0]*mb[2]
  m3=ma[0]*mb[1]-ma[1]*mb[0]
  return [m1,m2,m3]

#vector mode
def vcm(v):
  a=0
  for i in range(3):
    a+=(v[i])**2
  return(sqrt(a))

#vector scale
def vsc(v,s):
  vs=[v[0]*s,v[1]*s,v[2]*s]
  return vs

#distance of two points
def dist(a,b):
  d2=0
  for i in range(3):
    d2+=(a[i]-b[i])**2
  return(sqrt(d2))

#vector plus
def vplus(a,b):
  c=[]
  for i in range(3):
    c.append(a[i]+b[i])
  return c

#vector normalization
def normv(a):
  c=[]
  dv=1/vcm(a)
  for i in range(3):
    c.append(a[i]*dv)
  return(c)

#average of a set of vectors
def avec(A):
  n=len(A)
  s=[0,0,0]
  for a in A:
    for i in range(3):
      s[i]+=a[i]
  for i in range(3):
    s[i]=s[i]/n
  return(s)

#found projection of points on line
def projt(V,M,P):
  vp=vpts(M,P)
  v2=multp(V,V)
  bp=multp(V,vp)
  t=bp/v2
  x=[]
  for i in range(3):
    x.append(M[i]+t*V[i])
  return(x)

#vector connect two points
def vpts(a,b):
  vs=[b[0]-a[0],b[1]-a[1],b[2]-a[2]]
  return vs

#vector point product
def multp(ma,mb):
  return ma[0]*mb[0]+ma[1]*mb[1]+ma[2]*mb[2]

#vector mix product
def mmv(m,v):
  a1=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2]
  a2=m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2]
  a3=m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]
  return [a1,a2,a3]

#longitude and lattitude
def lltude(M):
  ll=[]
  lld=[]
  for a in M:
    xy=sqrt(a[0]**2+a[1]**2)
    phi=atan2(a[2],xy)
    theta=atan2(a[1],a[0])
    ll.append([theta,phi])
    thetad=theta*180/pi
    phid=phi*180/pi
    lld.append([thetad,phid])
  return(ll)

#distance matrix of SSE
def distm(crd):
#  for a in crd:
#    print a
  cas=len(crd)
#  print cas
  dis=[]
  i=0
  for i in range(cas):
#    print i
    dr=[]
    for k in range(cas):
#      print i,k
      x=0
      for j in range(3):
        t=crd[k][j]-crd[i][j]
        x+=t*t
      dr.append(sqrt(x))
    dis.append(dr)
#  i=0
#  for a in dis:
#    for a1 in a:
#      print str('%5.3f'%a1).zfill(6),
#    print
#    i+=1              
  return(dis)

def readin(f):
  f1=open(f,'r')
  fa=f1.readlines()
  f1.close()
  i=0
  for a in fa:
    if a[0]=='#':
      i+=1
    else:
      break
  fa=fa[i:]
  fa0=fa[:4]
#  print fa0
#  fd=fa0.split('\t')
  fd=[]
  for a in fa0:
    fd.append(a.strip())
  fd[1]=fa0[1][:-1]
  #print fd
  if fd[1].strip()=='':
    fd10=0
  elif fd[1][0]==' ':
    fd10=0
  else:
    fd10=int(fd[1][0])
#  if fd[2] in ['2','3','5']:
#    fd[2]=int(fd[2])
#  else:
#    print 'only 2,3,5 order of polynomial fit is acceptable; otherwise 2 is used'
#    fd[2]=2
  if fd10==2:
    if fd[3]=='':
      print ('no dssp program is given; use helix in pdb instead')   
      fd10=1
#   if fd10==3:
#     fexist=os.path.isfile(pathf+f)
#     if fexist==False:
#       print 'no ',name,' found; use helix in pdb instead'
#       fd10=1
#    elif fd[2]=='':
#      print 'no dssp program is given; use helix in pdb instead'   
#      fd10=1    
  plst=[]
  if fd[0]!='':
    if fd[0][-1]!='/':
      fd[0]+='/'
  if fd[2]!='':
    if fd[2][-1]!='/':
      fd[2]+='/'
  if fd[2]!='':    
    if not os.path.exists(fd[2]):
      os.mkdir(fd[2])
    
  for a in fa[4:]:
#    print a
    if a.strip()=='':
      continue
    b=a[:-1].split('\t')
    print (b)
    if b!=[]:
      plst.append(b)
  fd11=0
  if len(fd[1])>1:
    if fd[1][1]=='1':
      fd11=1
  if fd10==0:
    fd11=0
  fd[1]=[fd10,fd11]
#  print fd
  return([fd,plst])
