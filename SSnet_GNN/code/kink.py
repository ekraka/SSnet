#kink detection
from math import *

#irregular position scan based on rise
def kkr(rss):
	kk=[]
	n=len(rss)
	for i in xrange(n):
		rsi=rss[i]
		nr=len(rsi)
		kkh=[]
		if nr<7:
			kk.append(kkh)
			continue
		kav=sum(rsi)/nr
		for j in xrange(2,nr-2):
			r1=rsi[j]
			r0=abs(r1-rsi[j-1])
			r2=abs(r1-rsi[j+1])
			r10=r0/r1
			r20=r2/r1
			deltar=abs(r1/kav-1)
			if r1>2.:
				kkh.append([j,r1,deltar])
			elif deltar>0.67:
				kkh.append([j,r1,deltar])

# 			elif r10>0.75:
# #			elif r10>=4 or r10<=.25:
# 				kkh.append(j)
# 			elif r20>0.75:
# #			elif r20>=4 or r20<=.25:
# 				kkh.append(j)
		kk.append(kkh)
	return(kk)

def rish(ris,pdbh,nca):
#	print len(ris),nca
	rss=[]
	rn=[]
	for a in pdbh:
		rs=[]
		b0=a[0]-2
		b1=a[1]-2
		for i in xrange(b0,b1):
			rs.append(ris[i])
		rss.append(rs)
		rn.append(b0)
#	print 'rss',rss
	kd=kkr(rss)
#	print kd
	kkd=[]
	nh=len(kd)
	for i in xrange(nh):
		kkd1=[]
		if kd[i]==[]:
			kkd.append([])
			continue
		for a in kd[i]:
#			print rn[i]
#			print a
			kkd1.append([rn[i]+a[0]+1,a[1],a[2]])
		kkd.append(kkd1)
	return(kkd)


	