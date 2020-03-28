#helix shape classification

def hclass(data):
	rav,dltrr,nk=data[0],data[1],data[2]
	if rav >=100:
		tc='L'
	elif rav >=30:
		tc='LC'
	else:
		tc='HC'
	if dltrr<=0.15:
		tr='R'
	elif dltrr<0.35:
		tr='IR'
	else:
		tr='SI'
	if rav<=40 and dltrr>0.35:
		if nk>6:
			kk=1
		else:
			kk=0
	else:
		kk=0
	return([tc,tr,kk])

#def getkr():

def gethclass(axft):
#	kr=getkr()
	af=axft[1]
	akt=axft[2]
	ht=axft[0]
	#  hxt=axft[3]
	nh=len(ht)
	hct=[]
	for i in range(nh):
#	for a in kr:
		ki=akt[i][0]
		ti=akt[i][1]
		nk=len(ki)
		kmax=max(ki)
		kmin=min(ki)
		deltk=kmax-kmin
		rmax=1./kmin
		rmin=1./kmax
		deltr=rmax-rmin
		radii=[]
		for a in ki:
			radii.append(1./a)
		rav=sum(radii)/nk
		kav=sum(ki)/nk
#		dltrr=deltr/rmax
#		dltkk=deltk/kmax
		dltrr=deltr/rav
		dltkk=deltk/kav
		hca=hclass([rav,dltrr,nk])
		hpara=[rav,rmax,deltr,dltrr,kav,kmax,deltk,dltkk]
		hct.append([hca,hpara])
		
	return(hct)
	



