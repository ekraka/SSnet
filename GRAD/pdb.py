#Process pdb file; locate chain, model and coordinates from pdb file

from sys import *
import os
from io import StringIO as cStringIO
#import cStringIO

def pdbp2(pdb,chain,pathp,fh):

  f1=open(pathp+pdb+'.pdb','r')
  fa=f1.readlines()
  f1.close()
#determine chain id
  if chain=='':
    i=0
    for a in fa:
      la=a.split()
      if la[0]=='ATOM' and la[2]=='CA':
        if i==0:
          ch0=a[21]
          i+=1
        elif a[21]==ch0:
#          print ch0,la[4]
          i+=1
        else:
          i=0
        if i>10:
          chain=ch0
          break
#retrieve atoms of given chain
  fc=[]
  res=[]
  car=[]
  for a in fa:
    if (a[:3]=='TER' or a[:6]=='ENDMDL') and fc!=[]:
#      print 'end model'
#      print car,res
#      print fc
      break
    if len(a)<30:
      continue
    #print a
    if a[21]==chain and (a[:4]=='ATOM' or a[:6]=='HETATM'):
      if res!=[]:
        if a[22:26]==res[-1][22:26]:
#          print a
          res.append(a)
          car.append(a[13:15])
#          print car
        elif 'CA' in car:
#          print car
          for b in res:
#            print 'b',b
            fc.append(b)
          res,car=[a],[a[13:15]]
        else:
          res,car=[a],[a[13:15]]
      else:
        res.append(a)
        car.append(a[13:15])

  if 'CA' in car:
 #   print car
    for b in res:
      fc.append(b)
#  for a in fc:
#    print a
#  print 'fc',fc
  if fc==[]:
    return('')          
  fa=fc
      
   
#  fa=fh.getvalue()
#  print fa
#  fa=fb.readlines()
#  fh.seek(0)
#  fa=fh.readlines()
#  print fa
  fb=[]
  for a in fa:
    if a[:3]=='TER':
      break
    elif a[17:20]!='HOH':
      fb.append(a)
  nb=len(fb)
#  print 'nb',nb
  for i in range(nb-1,0,-1):
#    print fb[i].split()[0]
#    print fb[i]
    if fb[i][:4]!='ATOM':
      if fb[i][:6]!='HETATM':
        del fb[i]
      else:
        print ('HETATM',pdb, fb[i])
        break 
#  for a in fb:
#    print a 
  nb=len(fb)   
  fa,ca=[],[]
  for i in range(nb):
    a=fb[i]
#    la=a.split()
    if a[13:15]=='CA':
      ca.append([i,a[22:27],a[16]])
#  print ca
  nca=len(ca)
  cb,alt=[],[]
  for i in range(nca-1):
    if ca[i][1]==ca[i+1][1]:
      cb.append(ca[i])
      alt.append(ca[i][1])
#  alt = list(set(alt))
  ncb=len(cb)
#  print 'alt',alt
#  print 'cb',cb
  for i in range(nb):

    a=fb[i]
    ai=a[22:27]
    if ai in alt:
      j=alt.index(ai)
      if a[16]!=' ':
#        print a
        if a[16]==cb[j][2]:
          b=a[:16]+' '+a[17:]
          fa.append(b)
      else:
#        print a
        fa.append(a)
    else:
      fa.append(a)
  if chain==' ':
    name=pdb
  else:
    name=pdb+chain 
  nat=len(fa) 
  for i in range(nat):
    if fa[i][:6]=='HETATM':
#      print fa[i]
      fa[i]=fa[i].replace('HETATM','ATOM  ')
      
  f2=open(name+'p.pdb','w')
  for a in fa:
    f2.write(a)
  f2.close()
  return(chain)


def clean(pdb):
  os.system('rm '+pdb+'s.pdb')
  os.system('rm '+pdb+'p.pdb')


#process pdb file of NMR models
def pdbm(pdb,ch,mo,pathp,fh):

  f1=open(pathp+pdb+'.pdb','r')
  fa=f1.readlines()
  f1.close()
  nf=len(fa)
  p0,pt='',0
  for i in range(nf):
    a=fa[i]
    if a[:5]=='MODEL':
#      print a
      if a.split()[1]==mo:
#      if a[12:14]==mo:
        p0=i+1
#        j=i+1
#        print p0
        break
  if p0=='':
    if mo==1:
      p0==0
    else:
      return('')
  if p0!='':
    for i in range(p0+1,nf):
      a=fa[i]
      if a[:3]=='TER' or a[:6]=='ENDMDL' or a[:3]=='END':
        pt=i
#        print pt
        break
#  print p0,pt
  if pt==0 and mo==1:
    pt=nf-1
  if p0=='' and mo==1:
    p0==0
  if pt==0 or p0=='':
    return('')
  fb=fa[p0:pt]
  fc=[]
  res=[]
  car=[]
  for a in fb:
#    if a[21]==chain and (a[:4]=='ATOM' or a[:6]=='HETATM'):
    if a[:4]=='ATOM' or a[:6]=='HETATM':
      if res!=[]:
        if a[22:26]==res[-1][22:26]:
#          print a
          res.append(a)
          car.append(a[13:15])
#          print car
        elif 'CA' in car:
#          print car
          for b in res:
#            print 'b',b
            fc.append(b)
          res,car=[a],[a[13:15]]
        else:
          res,car=[a],[a[13:15]]
      else:
        res.append(a)
        car.append(a[13:15])
  if 'CA' in car:
 #   print car
    for b in res:
      fc.append(b)
#  for a in fc:
#    print a
#  print 'fc',fc
  if fc==[]:
    return('')          
  fa=fc
      
   
  fb=[]
  for a in fa:
    if a[:3]=='TER':
      break
    elif a[17:20]!='HOH':
      fb.append(a)
  nb=len(fb)
#  print 'nb',nb
  for i in range(nb-1,0,-1):
#    print fb[i].split()[0]
#    print fb[i]
    if fb[i][:4]!='ATOM':
#      print fb[i]
      del fb[i]
    else:
      break  
  nb=len(fb)   
  fa,ca=[],[]
  for i in range(nb):
    a=fb[i]
#    la=a.split()
    if a[13:15]=='CA':
      ca.append([i,a[22:27],a[16]])
#  print ca
  nca=len(ca)
  cb,alt=[],[]
  for i in range(nca-1):
    if ca[i][1]==ca[i+1][1]:
      cb.append(ca[i])
      alt.append(ca[i][1])
#  alt = list(set(alt))
  ncb=len(cb)
#  print 'alt',alt
#  print 'cb',cb
  for i in range(nb):

    a=fb[i]
    ai=a[22:27]
    if ai in alt:
      j=alt.index(ai)
      if a[16]!=' ':
#        print a
        if a[16]==cb[j][2]:
          b=a[:16]+' '+a[17:]
          fa.append(b)
      else:
#        print a
        fa.append(a)
    else:
      fa.append(a)
  name=pdb+mo  
  f2=open(name+'p.pdb','w')
  for a in fa:
    f2.write(a)
  f2.close()
  return(mo)

    
