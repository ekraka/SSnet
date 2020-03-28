#read axis trace point crd
#from sys import *

def raxt(ax,he):
#  print 'raxt ax',ax
  dx=[]
  for a in he:
    s0=a[0]
    s1=a[1]  
#    s2=s0-1
#    s3=s1
    s2=s0-2
    s3=s1-1
#    dx1=[]
    dx1=ax[s2:s3]
#    print a,dx
    dx.append(dx1)
  ndx=len(dx)
 
  return(dx)
