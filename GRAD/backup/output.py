#HAXIS out put
from figure import *

def breakins(seq,bki,bks):
  seq1=seq[:bki[0]]
  for i in range(len(bki)-1):
    bkinst=[]
    for k in range(bks[k]):
      bkinst.append(' ')
    seq1+=bkinst
    seq1+=seq[bki[i]:bki[i+1]+1]
  bkinst=[]
  for k in range(bks[-1]):
    bkinst.append(' ')
  seq1+=bkinst
  seq1+=seq[bki[-1]:]
  return(seq1)

def outputh(name,axft,ris,ax,vtnb,skts,pre0,pdbcrd,hlxc,kr,kq,axft3,axft5,hha,resi,sktasa,dual,o2,patho):
#  if dual==0 or (dual==1 and o2==1):
  hmethod={1:'pdb source',2:'DSSP',3:'user defined',0:'HAXIS'}
  ht=axft[0]
  hp=[]
  r3=range(3)
  for a in ht:
    h1=range(a[0]-1,a[1])
    hp+=h1
  n1=len(resi)
  seq=pdbcrd[2]
  crd=pdbcrd[1]
  if o2==1:
    f2=open(patho+name+'.haxis','w')
    f2.write('#PROGRAMM  HAXIS, v0.1.0a (30/03/13)\n\n')

  #input parameter
    f2.write('#Options:\n')
  #  f2.write('Path to pdb files:'+pre)
    option=['Path to pdb files: ','Helix definition option: ']
    opts=pre0[1][0]
    #for a in pre0[:2]:
    f2.write(option[0]+pre0[0]+'\n')
    if len(pre0)>1:
      f2.write(option[1]+str(opts)+' ('+hmethod[opts]+')\n')
    if dual==1:
      f2.write('Extra analysis by using HAXIS internal definition: On\n')
    f2.write('\n')


    
  #seqence and helix marked as H

  #  print 'seq',seq
    resi0=resi[0].strip()
    leads=[]



  #  print 'ht',ht 

  #  hs=leads
    hs=[]
    for i in range(n1):
      if i in hp:
        hs.append('H')
      else:
        hs.append(' ')
    ha=[]
    for a in hha:
      h1=range(a[0]-1,a[1])
      ha+=h1
  #  hsa=leads
    hsa=[]
    for i in range(n1):
      if i in ha:
        hsa.append('H')
      else:
        hsa.append(' ')        

    if resi0.isdigit():
      resi0d=int(resi0)
      lead=(resi0d-1)
      for i in range(lead):
        leads.append(' ')
      bks,bki=[],[]
      rid=resi0d
      for j in range(1,n1):
        rj=resi[j].strip()
        if rj.isdigit():
          rjd=int(rj)
          if rjd-rid==1:
            j+=1
            rid=rjd
          else:
            bks.append(rjd-rid-1)
            bki.append(j)
            j+=1
            rid=rjd
        else:
          j+=1
          rid+=1
  #    print 'bks',bks,bki


      if bki!=[]:
        seq1=breakins(seq,bki,bks)
        hs=breakins(hs,bki,bks)
        hsa=breakins(hsa,bki,bks)
        # seq1=seq[:bki[0]]
        # for i in xrange(len(bki)-1):
        #   bkinst=[]
        #   for k in xrange(bks[k]):
        #     bkinst.append(' ')
        #   seq1+=bkinst
        #   seq1+=seq[bki[i]:bki[i+1]+1]
   #     # seq1+=seq[bki[-1]:]
  #      print seq
  #      print seq1
      else:
        seq1=seq
    else:
      seq1=seq
  #  print leads
  #  print seq1
    seq2=leads+seq1
    hs=leads+hs
    hsa=leads+hsa
    n1=len(seq2)
  #  print n1
    f2.write('#Protein: '+name+' ('+str(n1)+' residues; source PDB)\n\n')
    lw=80
    n2=n1/lw
    n3=n1%lw
    f2.write('#Protein sequence and Helix definition\nline1: sequence\nline2: Helix defintion input\n')
    if opts!=0:
      f2.write('line3: Helix definition by HAXIS\n\n')
  #  f2.write(18*' '+'20'+18*' '+'40'+18*' '+'60'+18*' '+'80\n') 
                                          
    for i in range(int(n2)):
      numl,numv=[],[]
   #   f2.write(18*' ')
      numk=80*i
      for k in range(8):
   #     numi=80*i
        numk+=10
        numl.append(numk)
        numlk=10-len(str(numk))
        numv.append(numlk)
        f2.write(numlk*' '+str(numk))
      f2.write('\n')
  #    f2.write(18*' '+20+(i*80)+18*' '+'40'+18*' '+'60'+18*' '+'80\n') 
      for j in range(lw):
        f2.write(seq2[i*lw+j])
      f2.write('\n')
      for j in range(lw):
        f2.write(hs[i*lw+j])
      f2.write('\n')
      if opts!=0:
        for j in range(lw):
          f2.write(hsa[i*lw+j])
        f2.write('\n')
      f2.write('\n')
    n3c=n3/10
    numk=80*n2
    for k in range(int(n3c)):
      numk+=10
      numlk=10-len(str(numk))
      f2.write(numlk*' '+str(numk))
    f2.write('\n')
    for j in range(n3):
      f2.write(seq2[int(n2*lw+j)])
    f2.write('\n')
    for j in range(n3):
      f2.write(hs[n2*lw+j])
    f2.write('\n')
    if opts!=0:
      for j in range(n3):
        f2.write(hsa[n2*lw+j])
      f2.write('\n')
    f2.write('\n\n')
    n1=len(seq)
    
  #C alpha

    resi=pdbcrd[3]
  #  for a in crd:
  #    print a
    f2.write('#Residue number, Residue, Cartesian coordinates of anchor point C_alpha,\n\
     protein backbone length s, curvature kappa, torsion tau\n')
    arcl=[]
    arcs=0
    for a in skts[0]:
      arcs+=a
      arcl.append(arcs)
    skts[0]=arcl
    f2.write('  id   res        x         y         z         s         k         t\n')
  #               2	T	    -9.669    -0.447     4.998     0.000     0.000     0.000   
    for i in range(n1):
      f2.write(resi[i]+'\t'+seq[i]+'\t')
      for j in r3:
        f2.write('%10.3f'%crd[i][j])
      for j in r3:
        f2.write('%10.3f'%skts[i][j])
      f2.write('\n')
    f2.write('\n')
  #  for a in skts:
  #    print a
    f2.write('\n#Frenet unit vector T_i, N_i, and B_i (Coordinates with regard to the origin in the pdb coordinate system)\n')
    f2.write(' id                 T_i                              N_i                              B_i\n')
    f2.write('          x          y          z          x          y          z          x          y          z\n')  
  #           '       2.304    -1.658    -3.329     0.000     0.000     0.000     0.000    -0.000     0.000

  #  for a in vtnb:
  #    print a
  #  print len(crd),len(skts)
  #  for a in vtnb:
  #    print len(a)
    n2=len(vtnb[0])
    for i in range(n2):
      f2.write(resi[i])
      for j in r3:
        for k in r3:
          f2.write('%11.3e'%vtnb[j][i][k])
      f2.write('\n')
    f2.write('\n\n')  

  #Axis
    f2.write('#Local axis coordinates (projection of C alpha on axis)\n')
    f2.write('#For helix, they are real axis; otherwise, they are virtual or local axis based on the generalized axis\n\
   definition by treating strands, loops, turns etc. as fragments of helices and apply the same procedure.\n')
    f2.write('  id      x         y         z\n')
  #  print 'ax'
    j=0
    for a in ax:
  #    print a
      j+=1
      f2.write(resi[j])
      for i in r3:
        f2.write('%10.3f'%a[i])
      f2.write('\n')

  #Rise
    f2.write('\n\n#Parameter a_i (translation per residue along axis)\n')
    #for a in ris:
    nr=len(ris)
    for i in range(nr):
      f2.write('%5s'%resi[i+1]+'%10.3f'%ris[i]+' '+seq[i+1]+'\n')
    f2.write('\n\n')

  else:
    f2=open(name+'.haxis','a')
    f2.write('\n'+80*'*'+'\n\n')
    f2.write('#Extra analysis by using HAXIS internal helices definition\n')
    opts=0
  
#helices
  f2.write('#Helices input from '+hmethod[opts]+': residue(Start) - residue(End)\n')
  i=0
  for a in ht:
    i+=1
    f2.write('#%2d :'%i)
    f2.write(str(resi[a[0]-1])+'\t'+str(resi[a[1]-1])+'\n')

#helices
  if opts!=0:
    f2.write('#Helices identified by HAXIS: residue(Start) - residue(End)\n')
    i=0
    for a in hha:
      i+=1
      f2.write('#%2d :'%i)
      f2.write(str(resi[a[0]-1])+'\t'+str(resi[a[1]-1])+'\n')

  f2.write('#Currently, HAXIS look for segments showing regular geometry parameters:\n\
#rise parameter a, radius, arc length and curvature peak (of fitted axis curve), and\n\
#angle of neighboring local axis directions.\n\
#Where the regular region ends are the start/ends of helix.\n')
  
#Axis fitted & kt
  af=axft[1]
  akt=axft[2]
  ht=axft[0]
  hxt=axft[3]
  hfits=af
  hakts=akt
  hakts3=axft3[2]
  hakts5=axft5[2]
  hfits3=axft3[1]
  hfits5=axft5[1]
#  print len(ax),len(af),len(akt[0]),len(akt[0][0]),len(ax)
#  print 'af'
#  for a in af:
#    print a
#  print 'akt'
#  for a in akt:
#    print a
  nht=len(hxt)
  f2.write('\n\n#Anchor points of helix axis (projection of C_alpha position on helix axis)\n')
  for i in range(nht):
    f2.write('#%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+'\n')
    j=0
    f2.write('   id    x       y       z\n')
    for a in hxt[i]:
      f2.write('%5s'%resi[ht[i][0]-1+j])
      for b in a:
        f2.write('%8.3f'%b)
      f2.write('\n')
      j+=1
  forder=len(hfits[0][0])-1
  xl=['x','y','z']
  f2.write('\n#Polynominal fit of the helix axis: coefficients a_0 + a_1 x+a_2 x^2+...a_n x^n\n')
#  f2.write('\n#Polynominal fit of the helix axis: coefficients a0->an. 2nd order\n')
  f2.write('\n#2nd order: a0->a2\n')
  for i in range(nht):
    f2.write('#%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+'\n')
#    f2.write('		a0			a1				a2			a3\n')
    j=0
    for a in hfits[i]:
      f2.write(xl[j]+' ')
      for b in a:
        f2.write('%14.6e'%b)
      f2.write('\n')
      j+=1
#  f2.write('\n#Polynominal fit of the helix axis: coefficients a0->an. 3rd order\n')
  f2.write('#3rd order: a0->a3\n')
  for i in range(nht):
    f2.write('#%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+'\n')
#    f2.write('		a0			a1				a2			a3\n')
    j=0
    for a in hfits3[i]:
      f2.write(xl[j]+' ')
      for b in a:
        f2.write('%14.6e'%b)
      f2.write('\n')
      j+=1
#  f2.write('\n#Polynominal fit of the helix axis: coefficients a0->an. 5th order\n')
  f2.write('#5th order: a0->a5\n')
  for i in range(nht):
    f2.write('#%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+'\n')
    j=0
#    f2.write('		a0			a1				a2			a3\n')
    for a in hfits5[i]:
      f2.write(xl[j]+' ')
      for b in a:
        f2.write('%14.6e'%b)
      f2.write('\n')
      j+=1
  f2.write('\n#Curvature kappa and torsion tau of the helix axis presented by polynomial fits of order 2, 3, 5\n')
  for i in range(nht):
    f2.write('\nHelix # '+str(i+1)+': ')
    f2.write('%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+'\n\n')
    f2.write('                2nd order                 3rd order                  5th order\n')
    f2.write('           kappa        tau           kappa         tau            kappa        tau\n')
    haktsi=hakts[i]
    haktsi3=hakts3[i]
    haktsi5=hakts5[i]
    k=0
    for j in range(len(haktsi[0])):
      f2.write('%5s'%resi[ht[i][0]-1+k])
      f2.write('%14.6e'%haktsi[0][j]+'%14.6e'%haktsi[1][j])
      f2.write('%14.6e'%haktsi3[0][j]+'%14.6e'%haktsi3[1][j])
      f2.write('%14.6e'%haktsi5[0][j]+'%14.6e'%haktsi5[1][j]+'\n')
      k+=1
#    for a in hakts[i]:
#      for b in a:
#        f2.write('%14.6e'%b)
#      f2.write('\n')
    
#helix class
  # f2.write('\n#helix classification and global geometry parameters\n')
  # f2.write('#   Helix   Shap Reg.     Rav        Rmax        deltR      deltR/R       kav        kmax       deltk      deltk/k  Type\n')
  # for i in xrange(nht):
  #   f2.write('%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+' ')
  #   hct=hlxc[i][0]
  #   hpara=hlxc[i][1]
  #   f2.write('%4s'%hct[0]+'%4s'%hct[1]+'\t')
  #   for a in hpara:
  #     f2.write('%12.4e'%a)
  #   if ht[i][0]=='L':
  #   	kc='L'
  #   elif hct[2]==1:
  #   	kc='K'
  #   else:
  #   	kc='C'
  #   f2.write('  '+kc+'\n')

  f2.write('\n#Helix classification and global geometry parameters\n\n')
  f2.write('#Helix  (beginning and ending residue #) \n\
#Shape:   LC: weakly curved, HC: highly curved, L: linear\n\
#Regularity: R: regular, IR: irregular, SI: strongly irregular\n\
#C_av: average radius of curvature, C_max: maximum radius of curvature\n\
#delta_C: Variation of radius of curvature\n\
#delta_C/C_av\n\
#k_av: average curvature, k_max: maximum curvature\n\
#delta_k: Variation of radius of curvature\n\
#delta_k/k_av\n\
#Type: conventional type L: linear, C: curved, K: kinked\n')
  f2.write('\n#       Helix  Shape Reg.    C_av       C_max       Delta_C  Delta_C/C_av\
   k_av        k_max      Delta_k   Delta_k/k_av Type\n')
#  f2.write('      Helix\n')
  kk=[]
  for a in kq:
    if a!=[]:
      kk.append(a[-1]-1)
#  print kk,len(kr),nht
  for i in range(nht):
    f2.write('%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+' ')
    hct=hlxc[i][0]
    hpara=hlxc[i][1]
    f2.write('%4s'%hct[0]+'%4s'%hct[1]+'\t')
    for a in hpara:
      f2.write('%12.4e'%a)
    if ht[i][0]=='L':
      kc='L'
    elif hct[2]==1:
      if i in kk or kr[i]!=[]:
        kc='K'
      else:
        kc='C'
    else:
      kc='C'
    f2.write('   '+kc+'\n')

  # for i in xrange(nht):
  # f2.write('%4s'%hct[0]+'%4s'%hct[1]+'\t')
  # f2.write(
  # f2.write(`Helix  `)
  #   np=len(hpara)
  #   for i in xrange(np):
    
#kinks
  ck=0
  for a in kr:
    if a!=[]:
      ck+=1
  if ck!=0:
    f2.write('\n#Kinks (based on the irregularity of the rise parameter a_i)\n')
    f2.write('\n#       Helix     Kinks    a_i      delt_a/a_av\n')
    for i in range(nht):
      kri=kr[i]
      if kri==[]:
        continue
      f2.write('%2d: '%(i+1)+'%5s'%resi[ht[i][0]-1]+'%5s'%resi[ht[i][1]-1]+' : ')
      kri0=kri[0]
      f2.write('%5s'%resi[kri0[0]])
      for a in kri0[1:]:
        f2.write('%12.4e'%a)
      f2.write('\n')
      if len(kri)>1:

        for a in kri[1:]:
          f2.write(17*' '+'%5s'%resi[a[0]])
          for a1 in a[1:]:
            f2.write('%12.4e'%a1)
          f2.write('\n')

#kinkd based on k quartet
  hpk=kq
  ko=0
#  print 'outhpk',hpk
  if len(hpk)==[]:
    ko=1
#    return()
  else:
    i=0
    for a in hpk:
      if a!=[]:
        i+=1
    if i==0:
#      return()
      ko=1
      
#    f2=open(pdb+'.kpd2','w')
  if ko==0:
    f2.write('\n#Kinks (based on the axis curvature kappa)\n')
    f2.write('\n#pdb   kink ----Helix----  kinkregion  -------------- k quartet ----------- | ------------ rise -------------\n')
    for a in hpk:
      if a==[]:
        continue
#      print a[0]
#      print a[2]
      f2.write(name+'%5s'%resi[a[0][1]]+' %2d: '%a[-1]+'%5s'%resi[a[2][0]]+'%5s'%resi[a[2][-1]]+'%5s'%resi[a[0][0]]+'%5s'%resi[a[0][-1]])
#      for j in range(4):
#        f2.write('%5d'%a[3][j])
      f2.write('%10.4f'%a[1][1][0]+'%10.4f'%a[1][1][1]+'%10.4f'%a[1][1][2]+'%10.4f'%a[1][1][3])
      for j in range(len(a[3])):
        f2.write('%10.4f'%a[3][j])
#  f2.close()
      f2.write('\n')
  f2.write('\n#END\n')
  f2.close()
  
#write figure data
#  if dual==0 or (dual==1 and o2==1):
  if o2==1:
    name=name
  else:
    name=name+'2'
#   f2=open(name+'.dap','w')
#   f2.write('#Parameter a_i (translation per residue along axis)\n')
#   for i in xrange(1,n1-2):
#     f2.write('%5s'%resi[i]+'%10.4f'%ris[i-1]+' '+seq[i]+'\n')
#   f2.close()
  
  f2=open(patho+name+'.dca','w')
  f2.write('id	res	x	y	z	s	k	t	ax	ay	az\n')
  hpa=set(hp)-set((0,n1-1))
  for i in range(n1):
    ia=i-1
    f2.write(resi[i]+'\t'+seq[i]+'\t')
    for j in r3:
      f2.write('%10.3f'%crd[i][j])
#    sktai=[sktasa[i],skts[i][1],skts[i][2]]
    for j in r3:
      f2.write('%10.3f'%skts[i][j])
#      f2.write('%10.3f'%sktai[j])
    if i in hpa:
      for j in r3:
        f2.write('%10.3f'%ax[ia][j])
    f2.write('\n')    
  f2.close()
  
#write gp script
  wsgp(name,resi,ht,sktasa,o2,patho)
  return()
