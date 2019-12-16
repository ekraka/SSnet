#calculate arc length, k,t of spline points

from utils import *


def arc3(cc,x,ns):
#  ns=200
  ds=float(1./ns)
  sn=0
  rskt=[]
  for i in range(ns):
    d1,d2,d3=[],[],[]
    t=ds*i
    si=0
    for j in range(3):
      d3t=x[j]*2
      d2t=d3t*t+cc[j][2]*2
      d1t=x[j]*t*t+2*cc[j][2]*t+cc[j][1]
      si+=(d1t)**2
      d1.append(d1t)
      d2.append(d2t)
      d3.append(d3t)
    di=[d1,d2,d3]
    kap=curv(di[0],di[1])
    torp=tors(di[0],di[1],di[2])
    si=sqrt(si)*ds
    sn+=si
    rskt.append([sn,kap,torp])
  return (rskt)
  
def arc3a(cc,x,ns):
#  ns=200
  nt=ns*10
  ds=float(1./nt)
  sn=0
  rskt=[]
  for i in range(nt):
    d1,d2,d3=[],[],[]
    t=ds*i
    si=0
    for j in range(3):
      d3t=x[j]*2
      d1t=x[j]*t*t+2*cc[j][2]*t+cc[j][1]
      si+=(d1t)**2
    si=sqrt(si)*ds
    sn+=si
    if (i+1)%10==0:
      for j in range(3):
        d3t=x[j]*2
        d2t=d3t*t+cc[j][2]*2
        d1t=x[j]*t*t+2*cc[j][2]*t+cc[j][1]
        d1.append(d1t)
        d2.append(d2t)
        d3.append(d3t)
      di=[d1,d2,d3]
      kap=curv(di[0],di[1])
      torp=tors(di[0],di[1],di[2])
      rskt.append([sn,kap,torp])
    else:
      continue
  return (rskt)

def arcd(cc,x,ns):
#  ns=200
  nt=ns*10
  ds=float(1./nt)
  sn=0
  rskt=[]
  crds=[]
  for i in range(nt):
    d1,d2,d3=[],[],[]
    t=ds*i
    si=0
    for j in range(3):
      d3t=x[j]*2
      d1t=x[j]*t*t+2*cc[j][2]*t+cc[j][1]
      si+=(d1t)**2
    si=sqrt(si)*ds
    sn+=si
    if (i+1)%10==0:
      crd=[]
      for j in range(3):
        d3t=x[j]*2
        d2t=d3t*t+cc[j][2]*2
        d1t=x[j]*t*t+2*cc[j][2]*t+cc[j][1]
        d1.append(d1t)
        d2.append(d2t)
        d3.append(d3t)
        t2=t*t
        t3=t2*t
        crd.append(t3*x[j]/3+cc[j][2]*t2+cc[j][1]*t+cc[j][0])
      di=[d1,d2,d3]
      kap=curv(di[0],di[1])
      torp=tors(di[0],di[1],di[2])
      rskt.append([sn,kap,torp])
      crds.append(crd)
    else:
      continue
  return ([rskt,crds])
  
def sktd3(dat):
  skt=[]
  spl=len(dat[0][0])
  skts=[]
  ns=20
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
    sps,spk,spt=[],[],[]
    sktp=arc3b(cc,x0,ns)
    skt.append([s,k,tor])
    skts.append(sktp)
  return skts
  
def sktd3a(dat):
  skt=[]
  spl=len(dat[0][0])
  skts=[]
  ns=20
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
    sps,spk,spt=[],[],[]
    sktp=arc3a(cc,x0,ns)
    skt.append([s,k,tor])
    skts.append(sktp)
    
  sa=[0.]
  s0=0
  sr,sk,st=[0.],[0.],[0.]
  for a in skts:
#    s0=0
    s1=a[-1][0]
    for b in a:
      s=b[0]+s0
      sa.append(s)
      sr.append(b[0])
      sk.append(b[1])
      st.append(b[2])
    s0=s1+s0
  result=[sr,sa,sk,st]
  return result
  
def sktdc(dat):
  skt=[]
  spl=len(dat[0][0])
  skts=[]
  crds=[]
  ns=20
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
    sps,spk,spt=[],[],[]
    sktp=arcd(cc,x0,ns)
    skt.append([s,k,tor])
    skts.append(sktp[0])
    crds.append(sktp[1])
    
  sa=[0.]
  s0=0
  sr,sk,st=[0.],[0.],[0.]
  for a in skts:
#    s0=0
    s1=a[-1][0]
    for b in a:
      s=b[0]+s0
      sa.append(s)
      sr.append(b[0])
      sk.append(b[1])
      st.append(b[2])
    s0=s1+s0
  crd2=[]
  for a in crds:
    for b in a:
      crd2.append(b)
  result=[[sr,sa,sk,st],crd2]
  return result

