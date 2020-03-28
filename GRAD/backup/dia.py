from tool import *
def diam(crd,ax,p):
  nc=len(crd)
  na=len(ax)
#  print nc,na
  diamtr=[]
#  f2=open(p+'.dia','w')
  for i in range(na):
    distac=dist(ax[i],crd[i+1])
#    print distac
    diamtr.append(distac)
#    f2.write('%10.4f'%distac+'\n')
#  f2.close()
  
  return(diamtr)
  
    
