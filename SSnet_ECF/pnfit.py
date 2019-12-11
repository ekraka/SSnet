#curvature/torsion calculation of polynormial fitted line of axis
#2nd, 3rd and 5th order polynormial
from math import *
from utils import *

#5th order
def pfit(fth,t):

  nf=len(fth)
  o=nf-1
  L=t
  A=fth
#  for a in A:
#    print a

  #for each points
  D1,D2,D3=[],[],[]
  kapa,tora=[],[]
  crds=[]
  for r in xrange(L):
  #each point has xyz coordindates
    d1,d2,d3=[],[],[]
    crd=[]
    t=r+1
    for i in range(3):
  #    for j in range(3):
        ca=A[i]
#        print ca
        d1.append(ca[1]+ca[2]*t*2+ca[3]*t*t*3+ca[4]*t**3*4+ca[5]*t**4*5)
        d2.append(ca[2]*2+ca[3]*t*6+ca[4]*t**2*12+ca[5]*t**3*20)
        d3.append(ca[3]*6+ca[4]*t*12*2+ca[5]*t**2*20*3)
        crd.append(ca[0]+ca[1]*t+ca[2]*t**2+ca[3]*t**3+ca[4]*t**4+ca[5]*t**5)
    D1.append(d1)
    D2.append(d2)
    D3.append(d3)
 #   print crd
  #  print d1

    dcf=[]
    kap=curv(d1,d2)
    tor=tors(d1,d2,d3)
    kapa.append(kap)
    tora.append(tor)

#  print 'cur'
#  for a in kapa:
#    print a
#  print 'tor'
#  for a in tora:
#    print a
  return([kapa,tora])

#3rd order
def pfit3(fth,t):
  nf=len(fth)
  o=nf-1

  L=t
  A=fth
#  for a in A:
#    print a

  #for each points
  D1,D2,D3=[],[],[]
  kapa,tora=[],[]
  crds=[]
  for r in xrange(L):
  #each point has xyz coordindates
    d1,d2,d3=[],[],[]
    crd=[]
    t=r+1
    for i in range(3):
  #    for j in range(3):
        ca=A[i]
#        print ca
        d1.append(ca[1]+ca[2]*t*2+ca[3]*t*t*3)
        d2.append(ca[2]*2+ca[3]*t*6)
        d3.append(ca[3]*6)
        crd.append(ca[0]+ca[1]*t+ca[2]*t**2+ca[3]*t**3)
    D1.append(d1)
    D2.append(d2)
    D3.append(d3)
 #   print crd
  #  print d1

    dcf=[]
    kap=curv(d1,d2)
    tor=tors(d1,d2,d3)
    kapa.append(kap)
    tora.append(tor)

#  print 'cur'
#  for a in kapa:
#    print a
#  print 'tor'
#  for a in tora:
#    print a
  return([kapa,tora])

#2nd order
def pfit2(fth,t):

  nf=len(fth)
  o=nf-1

  L=t
  A=fth
#  for a in A:
#    print a

  #for each points
  D1,D2,D3=[],[],[]
  kapa,tora=[],[]
  crds=[]
  for r in xrange(L):
  #each point has xyz coordindates
    d1,d2,d3=[],[],[]
    crd=[]
    t=r+1
    for i in range(3):
  #    for j in range(3):
        ca=A[i]
#        print ca
        d1.append(ca[1]+ca[2]*t*2)
        d2.append(ca[2]*2)
        d3.append(0.)
        crd.append(ca[0]+ca[1]*t+ca[2]*t**2)
    D1.append(d1)
    D2.append(d2)
    D3.append(d3)
#    print crd
  #  print d1

    dcf=[]
    kap=curv(d1,d2)
    tor=tors(d1,d2,d3)
    kapa.append(kap)
    tora.append(tor)

#  print 'cur'
#  for a in kapa:
#    print a
#  print 'tor'
#  for a in tora:
#    print a
  return([kapa,tora])



