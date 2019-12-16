#helix kink detection from spline fit of axis traces
#curvature quartet

from math import *
from copy import *
from sys import *

def screen(ris,cut):
  he=[]
  n=len(ris)
  for i in xrange(n):
    if ris[i]<=cut:
      he.append(i)
  return(he)
  
def screenb(ris,cut):
  he=[]
  n=len(ris)
  for i in xrange(n):
    if ris[i]>=cut:
      he.append(i)
  return(he)
  
def screen2(ris,cutl,cuth):
  he=[]
  n=len(ris)
  for i in xrange(n):
    if ris[i]<=cuth and ris[i]>=cutl:
      he.append(i)
  return(he)

def segfp(he):
  herl=[]
  he1=[he[0]]
  nh=len
  for i in xrange(1,len(he)):
    if he[i]-he[i-1]==1:
      he1.append(he[i])
    elif len(he1)>2:
      herl.append(he1)
      he1=[he[i]]
    else:
      he1=[he[i]]
  if len(he1)>2:
    herl.append(he1)
  return(herl)


def hhrj(H,sak,ris,name,resi):
  pdb=name
  pkcut=1.0
  bkcut=0.4
  na=20
  he=[]
  for a in H:
    he1=range(a[0]-1,a[1])
    he.append(he1)
#   print 'helix'
#   for a in he:
#     print a
  ds=ris
#  ds.insert(0,0.)
#  ds.append(0.)
  nds=len(ds)
  #print nfds,nds,nfd
#  print nds

  kap=[]
  kac=[]
  kab=[sak[0]]
  n=len(sak)
  ka1=[sak[0]]
  for i in xrange(1,n-1):
    if i%20!=0:
      ka1.append(sak[i])
    else:
      kac.append(ka1)
      ka1=[sak[i]]
      ka1m=max(ka1)
      kap.append(ka1m)
      kab.append(sak[i])
  kab.append(sak[-1])
#  print len(kab),len(kap),len(kac)
#  print 'kab'
#  for a in kab:
#    print a
#  print 'kap'
#  for a in kap:
#    print a
#  for a in kac:
#    print len(a)
  kam=kap
  kap=kab
#  print 'he',he
#  kap.append(max(k1))
  nka=len(kap)
#  print kap
  ht=[]
  hpk=[]
  if he!=[]:
    for a in he:
      ca=deepcopy(a)

  #adjust entry point    
  #    for b in a:
  #      print b, ds[b]
      if ds[a[0]]<=2.:
        for i in range(1,3):
          id1=a[0]-i
          if id1<0:
            break
          if ds[id1]<2:
            a.insert(0,id1)
          else:
            break
      else:
        del a[0]
        for i in range(1,3):
          id2=a[0]+i
          if id2>=nds:
            break
          if ds[id2]>2:
            del a[0]
          else:
            break

  #adjust exit point
      if len(a)>=4:
        for i in range(2):
          if a[-1]>=nds:
            del a[-1]
          else:
            break  
        if ds[a[-1]]<=2.:
          for i in range(1,3):
            id1=a[-1]+i
            if id1<0 or id1>=nds:
              break
            if ds[id1]<2:
              a.append(id1)
            else:
              break
        else:
          del a[-1]
          for i in range(2):
#            print a
            id2=a[-1]
            if id2>=nds:
              break        
            if ds[id2]>2:
              del a[-1]
            else:
              break
      else:
        continue        
      
#      print a
      nba=len(a)
      if a[0]<ca[0]:
        d0=a[0]
      else:
        d0=ca[0]
      if a[-1]<ca[-1]:
        d1=ca[-1]
      else:
        d1=a[-1]
      nba=len(ca)
#       print ca
#       print a
#       for b in ca:
#         print resi[b],
#       print
#       for b in a:
#         print resi[b],
#       print

def hhs(H,sak,ris,name,resi,skts,diamt,axs):
  pdb=name
  pkcut=1.0
  bkcut=0.4
  rcut1=2.8
  rcut2=2.
  na=20
  
#rise
  he=[]
  hes=[]
  n=len(ris)
#  print n
#  heu=[]
#  for i in xrange(n):
#    if ris[i]<=rcut2:
#      he.append(i)
#    if ris[i]<rcut1:
#      heu.append(i)
  he=screen(ris,rcut2)
  heu=screen(ris,rcut1)
      

#  herl=[]
#  he1=[he[0]]
#  nh=len
#  for i in xrange(1,len(he)):
#    if he[i]-he[i-1]==1:
#      he1.append(he[i])
#    elif len(he1)>2:
#      herl.append(he1)
#      he1=[he[i]]
#    else:
#      he1=[he[i]]
#  if len(he1)>2:
#    herl.append(he1)

#  herh=[]
#  he1=[heu[0]]
#  nh=len
#  for i in xrange(1,len(heu)):
#    if heu[i]-heu[i-1]==1:
#      he1.append(heu[i])
#    elif len(he1)>2:
#      herh.append(he1)
#      he1=[heu[i]]
#    else:
#      he1=[heu[i]]
#  if len(he1)>2:
#    herh.append(he1)    
  herl=segfp(he)
  herh=segfp(heu) 
  
#arc length of axis spline
  scutl=2.
  scuth=3.
  arc=skts[0]
  harcl=screen(arc,scutl)
  harch=screen(arc,scuth)
  harcls=segfp(harcl)
  harchs=segfp(harch)
  
#arc length of axis spline
  kcut=1.
  cur=skts[1]
  hkc=screen(cur,kcut)
  hkcs=segfp(hkc)
  
#helix diameter
  dcutl=1.8
  dcuth=2.6
  hdm=screen2(diamt,dcutl,dcuth)
  hdms=segfp(hdm)

#vec angle
  ax=axs[0]
  axa=axs[1]
  axar=axs[2]
#  f2=open(name+'.axa','w')
  nax=len(axa)
#  for a in axa:
#  for i in xrange(nax):
#    f2.write('%10.4f'%axa[i]+'%10.4f'%axar[i]+'\n')
#  f2.close()
  acutl=-0.4
  acuth=0.
#  ham=screen2(axar,acutl,acuth)
  ham=screenb(axar,acutl)
  hams=segfp(ham)
  
# #ax direction angle degree
# #  ax=axs[0]
#   axb=axs[3]
#   axbr=axs[4]
# #  f2=open(name+'.axa','w')
#   nab=len(axb)
# #  for a in axa:
# #  for i in xrange(nax):
# #    f2.write('%10.4f'%axa[i]+'%10.4f'%axar[i]+'\n')
# #  f2.close()
#   bcutb=25.
#   bcutn=30.
# #  ham=screen2(axar,acutl,acuth)
#   hbm=screen(axb,bcutb)
#   hbms=segfp(hbm) 
  
#ax direction angle rad
  axb=axs[3]
  axbr=axs[4]
  ma=axs[5]
  nab=len(axb)
  bcutb=0.9
  bcutn=0.85
  if ma==0:
    bcut=bcutn
  elif ma==1:
    bcut=bcutb
#  ham=screen2(axar,acutl,acuth)
  hbm=screenb(axbr,bcut)
  hbms=segfp(hbm) 
      
#   print 'pdb'
#   for a in H:
#     print a
#   print 'haxis rl'
#   for a in herl:
#     print a
#   print 'haxis rh'
#   for a in herh:
#     print a
#   print 'haxis sl'
#   for a in harcls:
#     print a
#   print 'haxis sh'
#   for a in harchs:
#     print a
#   print 'haxis acur'
#   for a in hkcs:
#     print a
#   print 'haxis diam'
#   for a in hdms:
#     print a
#   print 'haxis vector angle'
#   for a in hams:
#     print a
#   print 'haxis ax2 angle'
#   for a in hbms:
#     print a
  
#   for a in [sak[0],ris,resi,skts[1],diamt,axs[2]]:
#     print len(a),
#   print
  n=len(resi)
  hit=[]
  for i in xrange(n):
    hit.append(0)
  for a in herl:
    for b in a:
      hit[b+2]+=1
    hit[b+3]+=1
  for a in herh:
    for b in a:
      hit[b+2]+=1
    hit[b+3]+=1
  for a in harcls:
    for b in a:
      hit[b+1]+=1
    hit[b+2]+=1
  for a in harchs:
    for b in a:
      hit[b+1]+=1
    hit[b+2]+=1
  for a in hdms:
    for b in a:
      hit[b+1]+=1
    hit[b+2]+=1
  for a in hams:
    for b in a:
      hit[b+1]+=1
    hit[b+2]+=1
  for a in hbms:
    for b in a:
      hit[b+1]+=1
    hit[b+2]+=1
    hit[b+3]+=1
#  for a in hit:
#  print hkcs
  hkcsa=[]
  for a in hkcs:
    hkcsa+=a
#  print 'hit',hit
  for i in xrange(n):
    if i-1 in hkcsa and hit[i]!=0:
      hit[i]+=1

#  print 'hit',hit
  his=screenb(hit,4)
#  print 'his',his
  hh=segfp(his)
#   print 'hh',hh
#   for a in hh:
#     print a
  hhe=[]
  for a in hh:
    hhe.append([a[0],a[-1]])
#  print hhe
#  print 'H',H
  return(hhe)
       

