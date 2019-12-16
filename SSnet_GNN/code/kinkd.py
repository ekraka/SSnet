#helix kink detection from spline fit of axis traces
#curvature quartet

from math import *
from copy import *
from sys import *

def kinkq(H,sak,ris,name,sktc):
  pdb=name
  pkcut=.8
  bkcut=0.25
  scut=2.5
  na=20
  nb=30
  nb2=80
  he=[]
  for a in H:
    he1=range(a[0]-1,a[1])
    he.append(he1)
#  for a in sktc:
#    print a
#  print 'helix'
#  for a in he:
#    print a
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
#  print he
  hi=0
  if he!=[]:
    for a in he:
      hi+=1
      if len(a)<7:
        continue
      ca=deepcopy(a)
# """
#   #adjust entry point    
#   #    for b in a:
#   #      print b, ds[b]
#       if ds[a[0]]<=2.:
#         for i in range(1,3):
#           id1=a[0]-i
#           if id1<0:
#             break
#           if ds[id1]<2:
#             a.insert(0,id1)
#           else:
#             break
#       else:
#         del a[0]
#         print a
#         for i in range(1,3):
#           id2=a[0]+i
#           if id2>=nds:
#             break
#           if ds[id2]>2:
#             del a[0]
#           else:
#             break
# 
#   #adjust exit point
#       if len(a)>=4:
#         for i in range(2):
#           if a[-1]>=nds:
#             del a[-1]
#           else:
#             break  
#         if ds[a[-1]]<=2.:
#           for i in range(1,3):
#             id1=a[-1]+i
#             if id1<0 or id1>=nds:
#               break
#             if ds[id1]<2:
#               a.append(id1)
#             else:
#               break
#         else:
#           del a[-1]
#           for i in range(2):
# #            print a
#             id2=a[-1]
#             if id2>=nds:
#               break        
#             if ds[id2]>2:
#               del a[-1]
#             else:
#               break
#       else:
#         continue        
#       
# #      print a
#       nba=len(a)
#       if a[0]<ca[0]:
#         d0=a[0]
#       else:
#         d0=ca[0]
#       if a[-1]<ca[-1]:
#         d1=ca[-1]
#       else:
#         d1=a[-1]
# """
      nba=len(ca)
#      print ca
  #    for i1 in ca:
  ##      print kap[i1-1],ds[i1-1]
  #    print
  #only look for kinks in helix minimum length 7
      if ca[-1]-ca[0]<=6:
        hpk.append([])
        continue
#      print 'ca',ca
      pk,bk,pki,bki=[],[],[],[]
      kmax,kmin,kmaxi,kmini=[],[],[],[]
#      knk=[]
      sseg=[]
      for l in xrange(3,nba-3):
        b=ca[l]-1
#        print kap[b]
        if kap[b]>=pkcut:
#          print 'peak',b,kap[b]
          pk.append(kap[b])
          pki.append(b)
#scan for k max/min. -1 to +4
          si=range(b-1,b+5)
          sseg+=si
      sseg=set(sseg)
      sseg=list(sseg&set(ca))
      sseg.sort()
#      print 'sseg',sseg
#merge scan sections
      sseg2=[]
      ns2=len(sseg)
      if ns2==0:
        continue
      ss1=[sseg[0]]
      for j in xrange(1,ns2):
        if sseg[j]-sseg[j-1]==1:
          ss1.append(sseg[j])
        else:
          sseg2.append(ss1)
          ss1=[sseg[j]]
      if ss1!=[]:
        sseg2.append(ss1)
#      print 'sseg2',sseg2
#scan
      for s1 in sseg2:
#        print s1
        maxs1,maxs1i=[],[]
        i0=s1[0]*na
        i1=(s1[-1]+1)*na
#        print i0,i1,sak[i0],sak[i1]
        for j in xrange(i0+1,i1):
          if sak[j]>sak[j-1] and sak[j]>sak[j+1] and sak[j]>pkcut:
            maxs1.append(sak[j])
            maxs1i.append(j)
#            print j,j/na
#        print maxs1,maxs1i
#search nwn peak pattern
        ns1=len(maxs1)
        if ns1 < 4:
          continue
        for j in xrange(ns1-3):
#          print 'search',j
          mj=maxs1i[j]
          mj1,mj2,mj3=maxs1i[j+1],maxs1i[j+2],maxs1i[j+3]
#          print mj,mj1,mj2,mj3
          if maxs1i[j+1]-maxs1i[j] <nb:
            if maxs1i[j+3]-maxs1i[j+2]<nb:
              maxsd1=maxs1i[j+2]-maxs1i[j+1]
              if maxsd1>na and maxsd1<nb2:
#                print 'qpeak find',mj/20,mj1/20,mj2/20,mj3/20
#                if mj3/20==ca[-1]:
#                  continue
                bk0,bk1,bk2='','',''
 #               bk1s=[]
 #               bk1is=[]
                bk1s={}
                for k in xrange(mj1,mj2):
#                  bk1s=[]
#                  bk1is=[]
                  if sak[k]<sak[k-1] and sak[k]<sak[k+1]:
                    bk1=sak[k]
#                    print bk1
                    bk1i=k
 #                   bk1s.append(bk1)
 #                   bk1is.append(bk1i)
                    bk1s[bk1]=bk1i
                if bk1s=={}:
                  continue
 #               print bk1s,bk1is
 #               print bk1s
 #               bk1m=min(bk1s)
 #               bk1mi=bk1is.index(bk1m)
 #               bk1m=sorted(bk1s.keys())[0]
                nbk1=len(bk1s)
                bk1m=sorted(bk1s.items(), key=lambda x: x[1])[0]
 #               print bk1m
 #               bk1mi=bk1s[bk1m]
                bk1=bk1m[0]
                bk1i=bk1m[1]
 #               bk1i=bk1is[bk1mi]
                for k in xrange(mj,mj1):
                  if sak[k]<sak[k-1] and sak[k]<sak[k+1]:
                    bk0=sak[k]
                    bk0i=k
                for k in xrange(mj2,mj3):
                  if sak[k]<sak[k-1] and sak[k]<sak[k+1]:
                    bk2=sak[k]
                    bk2i=k
#                print bk0,bk1,bk2,bk0i,bk1i,bk2i
                pq4=[maxs1i[j],maxs1i[j+1],maxs1i[j+2],maxs1i[j+3]]
#                print maxs1[j],maxs1[j+1],maxs1[j+2],maxs1[j+3],bk0,bk1,bk2,bk0i,bk1i,bk2i
                pqdat=[[maxs1i[j],maxs1i[j+1],maxs1i[j+2],maxs1i[j+3]],[maxs1[j],maxs1[j+1],maxs1[j+2],maxs1[j+3]],\
                [bk0i,bk1i,bk2i],[bk0,bk1,bk2]]
                if '' in [bk0,bk1,bk2]:
                  continue
                pq4a=[]
                pq4b=[]
                for pq1 in pq4:
                  pq4a.append(pq1/na)
                for pq1 in [bk0i,bk1i,bk2i]:
                  pq4a.append(pq1/na)
                  pq4b.append(pq1/na+1)
                  #                 if pq4b[-1]+1<=ca[-1]:
                pq4b.append(pq4b[-1]+1)
#                print pq4b
                if nbk1==1:
                  if bk1<bkcut and bk0>bkcut and bk2>bkcut:
                    pass
                  else:
                    continue
                elif nbk1>1:
                  if bk1<bkcut:
                  #                if bk1<bkcut and bk0>bkcut and bk2>bkcut:

                    s1,s2=sktc[0][pq4b[1]-1],sktc[0][pq4b[2]-1]
                  #                 print 's12',s1,s2, pq4b[1],pq4b[2]
                  #                 if s1<scut and s2<scut:
                    jsp=0
#                    print pq4b
                    for js in range(pq4b[0]+1,pq4b[-1]):
                      js1=sktc[0][js-1]
                  #                   print 'js',js1
                      if js1>scut:
                  	    jsp+=1
                    if jsp==0:
                      continue
#                  print pq4b
                  else:
                    continue
                else:
                  continue               
                if pq4b[-1] >=ca[-1]:
                  continue
                rd4=[]
                bl=range(pq4b[0],pq4b[-1])
#                  print 'bl',bl
                if len(bl)>3:
#                  for r1 in pq4b[:-1]:
                  for r1 in range(pq4b[0],pq4b[-1]):
                    rd1=ds[r1-1]
                    if rd1>2:
                      rd4.append(rd1)
                else:
                  for r1 in pq4b[:-1]:
                    rd1=ds[r1-1]
                    rd4.append(rd1)
#                  print rd4
#                  if rd4==[]:
#                    continue
#                  knk.append([pq4b,pqdat,ca,rd4])
#                print 'qpeak added',pq4b,ca
#                  if pq4b[-1] >=ca[-1]+1:
#                    continue
                hpk.append([pq4b,pqdat,ca,rd4,hi])
#                else:
#                  print 'not qpeak pattern'
#              else:
#                print 'not qpeak'
#            else:
#              print 'not qpeak'
#          else:
#            print 'not qpeak'
        kmax.append(maxs1)
        kmaxi.append(maxs1i)
      ht.append([a[0]+1,a[-1]+2])

#  for a in ht:
#    print a
#  print 'hpk',hpk
#  if len(hpk)==[]:
#    return()
#  else:
  i=0
  hpk2=[]
  # for a in hpk:
  #   print a
  # for a in H:
  #   print a
  for a in hpk:
    if a!=[]:
      i+=1
      hpk2.append(a)
  if i==0:
    return([])
  else:
#  if len(hpk)>0:
#    if hpk[0]==[]:
#      return([])
#    for a in hpk:
#      print a[0]
#      print a[2]    
    return(hpk2)

