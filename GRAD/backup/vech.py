#SSE vector analysis
#Vector directions, angles, dihedral angles, distances

from sys import *
from math import *
from tool import *
from rca import *
import os
from tors import *
#from output import *

def vecs(coef):
  patho=''
  patht=''
  pathpo=''
  n=len(coef[0][0])
#  print 'coef'
#  for a in coef:
#    print '#'
#    for b in a:
#      print b
#  print n,len(coef)
#  print n
  #T tangent,Y gama binormal, B beta normal
  T=[]
  for i in range(n):
    T.append([coef[2][0][i],coef[2][1][i],coef[2][2][i]])
#  print T
  T02=[]
  r3=range(3)
  for i in range(n):
    T01=[]
    for j in r3:
      t01=coef[0][j][i]*3+coef[1][j][i]*2+coef[2][j][i]
      T01.append(t01)
    T02.append(T01)
#   for i in xrange(n):
#     for j in r3:
#       print T[i][j],
#     for j in r3:
#       print T02[i][j],
#     print
  T.append(T02[-1])
#   print len(T)
#   print '#T'
#   for a in T:
#     print a
  
  Tn=[]
  for x in T:
    dv=0
    for i in range(3):
      dv+=x[i]*x[i]
    dv=sqrt(dv)
    td=[]
    for x1 in x:
      td.append(x1/dv)
    Tn.append(td)
#  print Tn

  #calculate normal and binormal vector
  D2=[]
  for i in range(n):
    d2=[]
    for j in range(3):
      d2.append(2*coef[1][j][i])
    D2.append(d2)
  D2a=[]
  for i in range(n):
    d2a=[]
    for j in r3:
      d2a1=coef[0][j][i]*6+coef[1][j][i]*2
      d2a.append(d2a1)
    D2a.append(d2a)
#  for x in D2:
#    print x
#  print '#Y'
  D2.append(D2a[-1])
#   print len(D2)
# 
#   for i in xrange(n):
#     for j in r3:
#       print D2[i][j],
#     for j in r3:
#       print D2a[i][j],
#     print
#   print len(D2)
  na=n+1
#   print 'D2'
#   for a in D2:
#     print a
  Y=[]
  for i in range(na):
    d12=multcros(T[i],D2[i])
    Y.append(d12)
#   for a in Y:
#     print a
#   print len(Y)
  
  Yn=[]
  for i in range(na):
    dv=vcm(Y[i])
#    print dv
    if dv!=0:
      Yn.append(vsc(Y[i],1/dv))
    else:
      Yn.append([0,0,0])
#   print 'Yn'
#   for a in Yn:
#     print a
#   print len(Yn)

#  print 'tangent vector'
#  printsf(T)
#  print 'normalized tangent vector'
#  printsf(Tn)
#   print 'Yn'
#   for a in Yn:
#     print a
#   print len(Yn)

  B=[]
  Bn=[]
  for i in range(na):
    d12=multcros(Y[i],T[i])
    B.append(d12)
  for i in range(na):
    db=multcros(Yn[i],Tn[i])
    Bn.append(db)
#   print 'normanl vector'
#   for a in B:
#     print a
#   print len(B)
# #  printsf(B)
#   print 'normalized normanl vector'
#   for a in Bn:
#     print a
#   print len(Bn)
  
#  printsf(Bn)
#  print 'binormal vector'
  #for a in Y:
  #  print a
#  printsf(Y)
#  print 'normalized binormal vector'
#  printsf(Yn)
  #TNB vector tuple for all residues(1~n-1) without two terminal residues.
  #rotation matrices A for all residues.
#  Tm=Tn[1:]
  A=[]
  for i in range(na):
    A.append([Tn[i],Bn[i],Yn[i]])
    
  nt=Tn
  nn=Bn
  nb=Yn
  M=[Tn,Bn,Yn]
  V=[T,B,Y]
  return(M)
#  return(V)

