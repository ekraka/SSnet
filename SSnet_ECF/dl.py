#generate axis traces from TNB vectors

from math import *
from sys import *
import numpy as np
#import cramer
from tool import *

def axis(vtnb,crd):
#axis method. ma=0: use normal; ma=1: use binormal and tangent
  ma=0
  
  VN=vtnb[1]
  nc=len(crd)
  nn=len(VN)
#  print nc,nn

#  print 'crd'
#  for a in crd:
#    print a
#   print 'normal v'
#   for a in VN:
#     print a

  def dot(p1, p2):
    return (p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2])

  def cross(p1, p2):
    return (p1[1]*p2[2] - p1[2]*p2[1], p1[2]*p2[0] - p1[0]*p2[2], p1[0]*p2[1] - p1[1]*p2[0])

  def pfactor(p, f):
    return (p[0]*f, p[1]*f, p[2]*f)

  def vplus(p1,p2):
    return(p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2])

  def vminus(p1,p2):
    return(p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2])

  def avv(p1,p2):
    return((p1[0]+p2[0])/2.,(p1[1]+p2[1])/2.,(p1[2]+p2[2])/2.)

  def mag(p):
    return sqrt(p[0]**2 + p[1]**2 + p[2]**2) 

  def normalize(p):
    m = mag(p)
    if m == 0:
      return(0.0, 0.0, 0.0)
    else:
      return(p[0]/m, p[1]/m, p[2]/m)

  def detme2(m):
    det=m[0][0]*m[1][1]-m[1][0]*m[0][1]
    return det

  def vdv(s1,s2):
    r3=range(3)
    d=[]
    for i in r3:
      d.append(s1[i]/s2[i])
    return(d)

  def linter(p1,v1,p2,v2):
    s1=cross(vminus(p2,p1),v2)
    s2=cross(v1,v2)
#    print s1,s2
    d=vdv(s1,s2)
    s3=sum(d)
    s3=s3/3
    #r3=range(3)
    #d=[]
    # s3=0
    # for i in r3:
    #   d1=s1[i]/s2[i]
    #   d.append(d1)
    #   s3+=d1
#    print s3/3
#    print d,s3
    l1=pfactor(v1,s3)
    l1=vplus(p1,l1)
    l2=vminus(l1,p2)
#    print l2,v2
    d2=vdv(l2,v2)
    s4=sum(d2)
    s4=s4/3
#    print l1,d2,s4
    return(l1)

  def lintera(p1,v1,p2,v2):
    m=[[v1[0],-v2[0],p2[0]-p1[0]],[v1[1],-v2[1],p2[1]-p1[1]]]
    c1=[v1[0],v1[1]]
    c2=[-v2[0],-v2[1]]
    c3=[p2[0]-p1[0],p2[1]-p1[1]]
    d=detme2([c1,c2])
    dx=detme2([c3,c2])
    dy=detme2([c1,c3])
    t1=dx/d
    t2=dy/d
    r3=range(3)
    cc1=[]
    for i in r3:
      cc1.append(v1[i]*t1+p1[i])
    cc2=[]
    for i in r3:
      cc2.append(v1[i]*t1+p1[i])   
#    print 'cc1',cc1
#    print 'cc2',cc2
    return(cc1)

  def cset(d1, d2):
    n=len(d1)
    crd=[d1[0]]
    for i in xrange(1,n):
        c1=d1[i]
        c2=d2[i-1]
        c3=[]
        for j in range(3):
            c=(c1[j]+c2[j])/2
            c3.append(c)
        crd.append(c3)
    crd.append(d2[-1])
    return(crd)

  #def lineq(m):
    

  ax=[]
  #axv average position for each segment (2 consecutive CA)
  axv=[]
  #axp average position for each residue (same residue for two consecutive segment)
  axp=[]
  axd=[]
  ax3=[]
  ax3v=[]
  #local ax direction angle of neighbor residue
  axar,axa=[],[]
  r3=range(3)
  for i in xrange(1,nn-2):
    u=normalize(VN[i])
    v=normalize(VN[i+1])
    #u=VN[0]
    #v=VN[1]
    ax1=cross(u,v)
#    print ax1
    axd.append(ax1)
    p=crd[i]
    q=crd[i+1]
#    print i,u,v,p,q
    a=dot(u,u)
    b=dot(u,v)
    axar.append(b)
    axa.append(acos(b)*180/pi)
    c=dot(v,v)
#    print a,b,c
    #w=cross(u,v)
    #wn=normalize(w)
    w=vminus(p,q)
#    print 'w',w,normalize(w)
    #w=wn
    d=dot(u,w)
    e=dot(v,w)
    D=a*c-b*b
#    print w,D
    sc=(b*e-c*d)/D
    tc=(a*e-b*d)/D
#    print sc,tc
    #dp=w+(sc*u)-(tc*v)
    #print sc,tc,dp
    P1=vplus(p,pfactor(u,sc))
    P2=vplus(q,pfactor(v,tc))
    ax.append([P1,P2])
    axv.append(avv(P1,P2))
#    print P1,P2
    v1=u
    v3=v
    v2=ax1
    L=[]
    for j in r3:
      L.append(v1[j])
      L.append(v2[j])
      L.append(-v3[j])
      L.append(q[j]-p[j])
    
#    A = np.array(L)
#    A.shape = (3,4)
#    result = cramer.solve(A)
#    if result:
#      x,y,z = result
#      print 'solution'
#      print 'x =', x
#      print 'y =', y
#      print 'z =', z, '\n'
#      cramer.check(A,x,y,z)
#     ax2i,ax4i=[],[]
#     ax3i=[]
#     for j in r3:
#       x2=p[j]+u[j]*result[0]
#       x4=q[j]+v[j]*result[2]
#       ax2i.append(x2)
#       ax4i.append(x4)
#       ax3i.append((x2+x4)/2.)
#       
#     ax3.append([ax2i,ax4i])
#     ax3v.append(ax3i)
  axp.append(ax[0][0])
  na=len(ax)
  for i in xrange(na-1):
    axp.append(avv(ax[i][1],ax[i+1][0]))
  axp.append(ax[-1][1])
#  return([axp,axa,axar])

  axtb=[]
  axn=[]
  A=vtnb
  axc=[]
  axtba=[]
  axtb1,axtb2=[],[]
  axb,axbr=[],[]
  for i in xrange(1,nn-2):
    B1=A[2][i]
    T1=A[0][i]
    B2=A[2][i+1]
    T2=A[0][i+1]
#    print i,B1,B2
    db=dist(B1,B2)
    dt=dist(T1,T2)
  #  db=dist(B1,B2)+dist(B1,Yn[i-1])
  #  dt=dist(T1,T2)+dist(T1,Tm[i-1])
    rbt=db/dt
    T1s=vsc(T1,rbt)
    AXtb=vplus(B1,T1s)
    AXtbn=normv(AXtb)
  #    a1=abs(ln[i][0]-ln[j][0])
  #    a2=abs(ln[i][0]+ln[j][0])
  #    a3=abs(a1-pi)
  #    if a1<pi/6 or a2<pi/6 or a3<pi/6:
  #      continue
  #    else:
#    cj+=1
#    na=multcros(A[1][i],A[1][i+1])
#    nan=normv(na)
#    axn.append(na)
    axtb.append(AXtbn)
    u=normalize(VN[i])
    v=normalize(VN[i+1])
    p=crd[i]
    q=crd[i+1]
    dtb1=pfactor(cross(AXtbn,cross(u,AXtbn)),1/dot(AXtbn,AXtbn))
    dtb2=pfactor(cross(AXtbn,cross(v,AXtbn)),1/dot(AXtbn,AXtbn))
#    print dtb1,dtb2
    qp=vminus(q,p)
    dqp=dot(qp,AXtbn)
    pjqp=vminus(q,pfactor(AXtbn,dqp))
#    caxtb=linter(p,dtb1,q,dtb2)
#    caxtb1=lintera(p,dtb1,q,dtb2)
    caxtb1=lintera(p,dtb1,pjqp,dtb2)
    pq=vminus(p,q)
    dpq=dot(pq,AXtbn)
    pjpq=vminus(p,pfactor(AXtbn,dpq))
    caxtb2=lintera(q,dtb2,pjpq,dtb1)
#    caxtb1=lintera(p,dtb1,q,dtb2)
#    print 'ax',caxtb1,caxtb2
    axtb1.append(caxtb1)
    axtb2.append(caxtb2)
#  print len(axtb1),len(axtb2)
#  return([axp,axa,axar])

  axtba=cset(axtb1,axtb2)
  n1=len(axtb)
  for i in xrange(n1-1):
    b=dot(axtb[i],axtb[i+1])
    axbr.append(b)
    axb.append(acos(b)*180/pi)
  axc,axcr=[],[]
  for i in xrange(n1-1):
    b=dot(axd[i],axd[i+1])
    axcr.append(b)
    axc.append(acos(b)*180/pi)  
#  print len(axa),len(axb),len(axc)
  n1=len(axa)
  n2=len(axb)
  n3=min(n1,n2)
#  for i in xrange(n3):
#    print axar[i],axbr[i],axcr[i]
  
#    print AXtbn,nan
#  print 'direction'
#  for a in axd:
#    print a

#  print 'axis direction 3 lines'
#  f2.write('axis direction 3 lines\n')
#  for a in ax3:
#    for b in a:
#      for c in b:
#        print c,
#        f2.write('%8.3f'%c)
#    print
#    f2.write('\n')

#  print 'axis direction 3 lines average'
#  f2.write('axis direction 3 lines average\n')
#  for a in ax3v:
#    for b in a:
#      print b,
#      f2.write('%8.3f'%b)
#    print
#    f2.write('\n')
#  f2.close()
#  data=[ax,axv,axp,axd,ax3,ax3v]
#   print 'axd'
#   for a in axd:
#    print a
#  n2=len(axd)
  # print '#########'
#  print len(axp),len(axtba)
#  n=len(axtba)
#  for i in xrange(n):
#    print mag(vminus(axtba[i],axp[i+1]))
  #   print axd[i],axtb[i]
#    print axtba[i],axp[i+1],mag(vminus(axtba[i],axp[i+1]))
#    for a in axtba[i]:
#    for a in axp[i+1]:
#    for a in ax[i][0]:
#     for a in axp[i]:
#       print a,
#     for a in axtba[i]:
#       print a,
#     print
  if ma==0:
    data=[axp,axa,axar,axc,axcr,ma]
  elif ma==1:
    data=[axtba,axa,axar,axb,axbr,ma]
  return(data)   
