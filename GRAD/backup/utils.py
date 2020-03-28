#basic functions module
from math import *
from copy import *
def bcs(b,c):
  return b+c
  
def bcsum(cc):
  bc=[]
  for i in range(3):
    bc.append(bcs(cc[i][1],cc[i][2]))
  return bc 
  
def multcros(ma,mb):
  m1=ma[1]*mb[2]-mb[1]*ma[2]
  m2=ma[2]*mb[0]-ma[0]*mb[2]
  m3=ma[0]*mb[1]-ma[1]*mb[0]
  return [m1,m2,m3]
  
def multp(ma,mb):
  return ma[0]*mb[0]+ma[1]*mb[1]+ma[2]*mb[2]
  
def multmix(ma,mb,mc):
  mbc=multcros(mb,mc)
  return multp(ma,mbc)

def detm(m):
  det=m[0][0]*m[1][1]*m[2][2]+m[0][1]*m[1][2]*m[2][0]+m[0][2]*m[1][0]*m[2][1]-m[2][0]*m[1][1]*m[0][2]-m[2][1]*m[1][2]*m[0][0]-m[2][2]*m[1][0]*m[0][1]
  return det

def invmd(m):
  det=detm(m)
  a11=(m[1][1]*m[2][2]-m[2][1]*m[1][2])/det
  a12=(m[0][2]*m[2][1]-m[2][2]*m[0][1])/det
  a13=(m[0][1]*m[1][2]-m[1][1]*m[0][2])/det
  a21=(m[1][2]*m[2][0]-m[2][2]*m[1][0])/det
  a22=(m[0][0]*m[2][2]-m[2][0]*m[0][2])/det
  a23=(m[0][2]*m[1][0]-m[1][2]*m[0][0])/det
  a31=(m[1][0]*m[2][1]-m[2][0]*m[1][1])/det
  a32=(m[0][1]*m[2][0]-m[2][1]*m[0][0])/det
  a33=(m[0][0]*m[1][1]-m[1][0]*m[0][1])/det
  return[[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]
  
def invm(m):
  a11=m[1][1]*m[2][2]-m[2][1]*m[1][2]
  a12=m[0][2]*m[2][1]-m[2][2]*m[0][1]
  a13=m[0][1]*m[1][2]-m[1][1]*m[0][2]
  a21=m[1][2]*m[2][0]-m[2][2]*m[1][0]
  a22=m[0][0]*m[2][2]-m[2][0]*m[0][2]
  a23=m[0][2]*m[1][0]-m[1][2]*m[0][0]
  a31=m[1][0]*m[2][1]-m[2][0]*m[1][1]
  a32=m[0][1]*m[2][0]-m[2][1]*m[0][0]
  a33=m[0][0]*m[1][1]-m[1][0]*m[0][1]
  return[[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]]

def mmv(m,v):
  a1=m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2]
  a2=m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2]
  a3=m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2]
  return [a1,a2,a3]

def drv(cc,x):
  d1=[]
  d2=[]
  d3=[]
  for i in range(3):
    d1.append(x[i]+2*cc[i][2]+cc[i][1])
    d2.append(2*x[i]+2*cc[i][2])
    d3.append(2*x[i])
  return([d1,d2,d3])

def arc(cc,x):
  ns=200
  ds=float(1./ns)
  sn=0
  for i in range(ns):
    t=ds*i
    si=0
    for j in range(3):
      si+=(x[j]*t*t+2*cc[j][2]*t+cc[j][1])**2
    si=sqrt(si)*ds
    sn+=si
  return (sn)
  
def arcj(cc,x):
  ns=200
  ds=float(1./ns)
  sn=0
  pa=[0.,0.,0.]
  for i in range(ns):
    t=ds*i
    si=0.
    ri=[]
    for j in range(3):
      ri.append(x[j]*t*t+2*cc[j][2]*t+cc[j][1])
      si+=ri[j]*ri[j]
#    print 'ri', ri
    sir=sqrt(si)
    sit=sir*ds
    sn+=sit
    for k in range(3):
      pa[k]+=(t*t*ds*ri[k])/sir
#    print 'pa', pa
  return ([sn,pa[0],pa[1],pa[2]])  
"""
def curv(d1,d2):
  d12=(d1[0])**2+(d1[1])**2+(d1[2])**2
  d22=(d2[0])**2+(d2[1])**2+(d2[2])**2
  cur=d22/(d12**2)
  return(cur)
"""
def curv(d1,d2):
  d12=(d1[0])**2+(d1[1])**2+(d1[2])**2
  d22=(d2[0])**2+(d2[1])**2+(d2[2])**2
  cur=(d22/(d12**2))**.5
  return(cur)

def tors(d1,d2,d3):
  tt=multmix(d1,d2,d3)
  tb=multp(d1,d1)*multp(d2,d2)
  if tb!=0:
    tor=(tt/tb)
  else:
    tor=0.
  return(tor)

def jacobfd(f1,x0,dx,skt,cc):
  aax=[]
  ax=[]
  for i in range(3):
    for j in range(3):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  for i in range(3):
    j1.append((arc(cc,aax[i])-skt[0]-f1[0])/dx)
    dj=drv(cc,aax[i])
    j2.append((curv(dj[0],dj[1])-skt[1]-f1[1])/dx)
    j3.append((tors(dj[0],dj[1],dj[2])-skt[2]-f1[2])/dx)
  J=[j1,j2,j3]
  return(J)
  
def arc0(xin0):
  ns=1000
  ds=float(1./ns)
  sn=0
  for i in range(ns):
    t=ds*i
    si=sqrt((xin0[0]*t*t+xin0[1])**2+(xin0[2]*t*t)**2)*ds
    sn+=si
  return (sn)

def curv0(x):
  cur=(4*x[0]**2+4*x[2]**2)/((x[0]+x[1])**2+x[2]**2)**2
  return(cur**.5)

def dst(x):
  dist=sqrt((x[0]/3+x[1])**2+(x[2]/3)**2)
  return dist
  
def jacobfd0(f1,x0,dx,skt,cbond):
  aax=[]
  ax=[]
  for i in range(3):
    for j in range(3):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  for i in range(3):
    j1.append((arc0(aax[i])-skt[0]-f1[0])/dx)
    j2.append((curv0(aax[i])-skt[1]-f1[1])/dx)
    j3.append((dst(aax[i])-cbond-f1[2])/dx)
  J=[j1,j2,j3]
  return(J)
  
def coefbd(c0):
  c1=[]
  for i in range(3):
    d=c0[i][0]+c0[i][1]+c0[i][2]+c0[i][3]
    c=c0[i][3]*3+c0[i][2]*2+c0[i][1]
    b=c0[i][3]*3+c0[i][2]
    c1.append([d,c,b])
  return c1
  
def repeatcb(xsol,x0,ct1,ct2,n,ito,st,cc,sovab):
    bs02=0
    lx=len(xsol)
    if lx>0:
#       ixc=0
#       bs02=0
        for xr in range(lx):
            bs01=0
            for ixr in range(3):
                if abs(xsol[xr][ixr]-x0[ixr])<1.e-5:
                    bs01+=1
            if bs01==3:
                ct1[xr][4]+=1
                axr=0
#                ax0=0
                for i in range(3):
                    axr*=abs(xsol[xr][i])/abs(x0[i])
                if axr>1:
                    xsol[xr]=deepcopy(x0)
                ct2.append([0,st,n,ito])
                bs02=0
                break
            else:
                bs02+=1
#         ixc+=1
        if bs02>1e-10:
            xsol.append(x0)
            ct1.append(['CA',0,st,n,1])
            ct2.append([0,st,n,ito])
    else:
            xsol.append(x0)
            sovab.append(cc)
            ct1.append(['CA',0, st,n,1])
            ct2.append([0,st,n,ito])
    return bs02

def jacobfd02(f1,x0,dx,skt,cbond):
  aax=[]
  ax=[]
  for i in range(3):
    for j in range(3):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  for i in range(3):
    j1.append((arc02(aax[i])-skt[0]-f1[0])/dx)
    j2.append((curv02(aax[i])-skt[1]-f1[1])/dx)
    j3.append((dst02(aax[i])-cbond-f1[2])/dx)
  J=[j1,j2,j3]
  return(J)
  
def arc02(xin0):
  ns=1000
  ds=float(1./ns)
  sn=0
  for i in range(ns):
    t=ds*i
    si=sqrt((xin0[0]*t*t+xin0[1])**2+(xin0[2]*t*t-xin0[2]/3)**2)*ds
    sn+=si
  return (sn)

def curv02(x):
  cur=(4*x[0]**2+4*x[2]**2)/((x[0]+x[1])**2+(2*x[2]/3)**2)**2
  return(cur**.5)

def dst02(x0):
  dist=x0[0]/3+x0[1]
  return dist

def jacobfd2d(f1,x0,dx,skt,cc):
  aax=[]
  ax=[]
  for i in range(2):
    for j in range(2):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  for i in range(2):
    j1.append((arc(cc,aax[i])-skt[0]-f1[0])/dx)
    dj=drv(cc,aax[i])
    j2.append((curv(dj[0],dj[1])-skt[1]-f1[1])/dx)
  J=[j1,j2]
  return(J)
  
def detm2(m):
  det=m[0][0]*m[1][1]-m[1][0]*m[0][1]
  return det
  
def invm2d(m):
  det=detm2(m)
  a11=m[1][1]/det
  a12=-m[0][1]/det
  a21=-m[1][0]/det
  a22=m[0][0]/det
  return[[a11,a12],[a21,a22]]

def jacobfd5(f1,x0,dx,skt,cc):
  aax=[]
  ax=[]
  for i in range(5):
    for j in range(5):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  j4=[]
  j5=[]
  for i in range(5):
    j1.append((arc(cc,aax[i])-skt[0]-f1[0])/dx)
    j2.append((curv(dj[0],dj[1])-skt[1]-f1[1])/dx)
    j3.append((arc(cc,aax[i])-skt[0]-f1[0])/dx)
    dj=drv(cc,aax[i])
    j4.append((curv(dj[0],dj[1])-skt[1]-f1[1])/dx)
    j5.append((tors(dj[0],dj[1],dj[2])-skt[2]-f1[2])/dx)
  J=[j1,j2,j3,j4,j5]
  return(J)
  
def sktd(dat):
  skt=[]
  spl=len(dat[0])
  for i0 in range(0, spl):
    cx=dat[0][i0]
    cy=dat[1][i0]
    cz=dat[2][i0]
    cx.insert(0,0)
    cy.insert(0,0)
    cz.insert(0,0)
    cc=[cx,cy,cz]
    x0=[3*cx[3],3*cy[3],3*cz[3]]
    bc=bcsum(cc)
    di=drv(cc,x0)
    s=arc(cc,x0)
    k=curv(di[0],di[1])
    if i0<spl-1:
      tor=tors(di[0],di[1],di[2])
    else:
      tor=0
    skt.append([s,k,tor])
  return skt

def jk(d1,d2,d12,d22):
  j2=[]
  for i in range(3):
#    j2.append((d12*d12*d2[i]*4-d22*d12*d1[i]*4)/(d12*d12*d12*d12))
    j2.append((d12*d2[i]*4-d22*d1[i]*4)/d12/d12/d12)
  return(j2)
  
def jacob(f1,x0,dx,skt,cc):
  aax=[]
  ax=[]
  t2=t*t
  t3=t2*t
  for i in range(3):
    for j in range(3):
      ax.append(x0[j])
    ax[i]+=dx
    aax.append(ax)
    ax=[]
  j1=[]
  j2=[]
  j3=[]
  for i in range(3):
    j1.append((arc(cc,aax[i])-skt[0]-f1[0])/dx)
    dj=drv(cc,aax[i])
    j2.append((curv(dj[0],dj[1])-skt[1]-f1[1])/dx)
    j3.append((tors(dj[0],dj[1],dj[2])-skt[2]-f1[2])/dx)
  J=[j1,j2,j3]
  return(J)
  
def readcrd(fname):
  f1=open(fname,'r')
  i=0
  pdb=[]
  crd=[]
  for line in f1:
    s=line.split()
    if s[0]=='TER':
      break
    else:
      if len(s)>7:
        if s[2]=='CA':
          if s[0]=='ATOM':
            pdb.append([float(s[6]),float(s[7]),float(s[8])])
  return(pdb)


def splcoef(x):
  n=len(x)
#  print n
  c=[0.]
  b=[4.]
  a=[1.]
  r=[0.]
  for i in range(1,n-2):
    c.append(1.)
    a.append(1.)
  for i in range(1,n-1):
    b.append(4.)
    r.append(3*(x[i+1]-2*x[i]+x[i-1]))
  s=b[0]
  t=[0.]
  u=[]
  n=n-1
  u.append(r[0]/s)
  for j in range(1,n):
    t.append(c[j-1]/s)
    s=b[j]-a[j-1]*t[j]
    u.append((r[j]-a[j-1]*u[j-1])/s)
  for j in range(n-2,-1,-1):
    u[j]-=t[j+1]*u[j+1]
  return(u)
  
def sktd2(dat):
  skt=[]
  spl=len(dat[0][0])
  for i0 in range(0, spl):
    cx=[dat[0][0][i0],dat[1][0][i0],dat[2][0][i0],dat[3][0][i0]]
    cy=[dat[0][1][i0],dat[1][1][i0],dat[2][1][i0],dat[3][1][i0]]
    cz=[dat[0][2][i0],dat[1][2][i0],dat[2][2][i0],dat[3][2][i0]]
#    cx.insert(0,0)
#    cy.insert(0,0)
#    cz.insert(0,0)
    cc=[cx,cy,cz]
    x0=[3*cx[3],3*cy[3],3*cz[3]]
    bc=bcsum(cc)
    di=drv(cc,x0)
    s=arc(cc,x0)
    k=curv(di[0],di[1])
    if i0<spl-1:
      tor=tors(di[0],di[1],di[2])
    else:
      tor=0
    skt.append([s,k,tor])
  return skt

  
