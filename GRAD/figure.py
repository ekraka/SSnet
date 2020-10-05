def gsect(p,la,lb,f2,wd,patho):

	f2.write('set multiplot\n\
set size '+str(wd)+',0.475\n')


	f2.write('set origin 0.0,0.5          # kappa(s)\n'
	+'set xrange ['+str(la)+':'+str(lb)+']\n'
	+'set yrange [0:4]\n')


	f2.write('plot "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0) \\\n\
	  w l lt 1 lw 1 title "{/Symbol k}", \\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0):10 \\\n\
	  w labels rotate left offset 0,1 notitle, \\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0):11 \\\n\
	  w labels rotate left offset 0,2 notitle, \\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $8 : 1/0) \\\n\
	  pointtype 7 pointsize 0.5 notitle\n\n')


	f2.write('set origin 0.0,0.0          # tau(s)\n\
set yrange [-1:1]\n')
	f2.write('plot "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $4 : 1/0) \\\n\
	  w l lt 1 lw 1 title "{/Symbol t}", \\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $4 : 1/0):10 \\\n\
	  w labels rotate left offset 0,1 notitle, \\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $4 : 1/0):11 \\\n\
	  w labels rotate left offset 0,2 notitle ,\\\n\
	  "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $9 : 1/0) \\\n\
	  pointtype 7 pointsize 0.5 notitle\n\n')

	f2.write('unset multiplot\n\n')

def gseck(p,la,lb,f2,wd,og):
  f2.write('set origin 0.0,'+og+'          # kappa(s)\n'
  +'set xrange ['+str(la)+':'+str(lb)+']\n'
  +'set yrange [0:4]\n'
  +'set size '+str(wd)+',0.475\n')


  f2.write('plot "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0) \\\n\
    w l lt 1 lw 1 title "{/Symbol k}", \\\n\
    "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0):10 \\\n\
    w labels rotate left offset 0,1 notitle, \\\n\
    "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $3 : 1/0):11 \\\n\
    w labels rotate left offset 0,2 notitle, \\\n\
    "'+p+'.axsp" i 0 u 2:($2 >= '+str(la)+' && $2 < '+str(lb)+' ? $8 : 1/0) \\\n\
    pointtype 7 pointsize 0.5 notitle\n\n')
	
def rsec(p,b0,b1,f2,wd):
    f2.write('set size '+str(wd)+',0.475\n')
    f2.write('set xrange ['+str(b0)+':'+str(b1)+']\n')
    f2.write('set yrange [0:10]\n')
    f2.write('plot "'+p+'.haxis" i 5 u 1:($1 >='+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0) \\\n\
	  w l lt 1 lw 1 title "{a}", \\\n\
	  "'+p+'.haxis" i 5 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0):3 \\\n\
	  w labels rotate left offset 0,1 notitle, \\\n\
	  "'+p+'.haxis" i 5 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0):1 \\\n\
	  w labels rotate left offset 0,2 notitle, \\\n\
	  "'+p+'.haxis" i 5 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0) \\\n\
	  pointtype 7 pointsize 0.5 notitle\n\n')

def r3d(p,b0,b1,f2):
    f2.write('unset multiplot\n\n\
  ### 3D structure ###\n\
  set origin 0,0; set size 1,1\n\
  set autoscale\n\
  set format "%2.0f"\n\
  set xtics autofreq; set ytics autofreq; set ztics autofreq\n\
  set xlabel "x"; set ylabel "y"; set zlabel "z"\n\
  set border 4095\n\
splot "'+p+'.dca" every ::'+str(b0)+'::'+str(b1)+' i 0 u 3:4:5 w l lt 1 title "r(s)", \\\n\
        "'+p+'.dca" every ::'+str(b0)+'::'+str(b1)+' i 0 u 3:4:5:2 w labels notitle,\\\n\
 "'+p+'.dca" every ::'+str(b0)+'::'+str(b1)+' i 0 u 9:10:11 \\\n\
    pointtype 7 pointsize 0.5 lc rgb "red" notitle\n')

def wsgp(p,resi,hp,sa,o2,patho):
  p1=p
  if o2==2:
    p1=p[:-1]
  n=int(resi[-1])
  nr=n/50+1
  nr2=nr/2
  nra=nr%2
  nrb=n%50
  we=float(nrb)/50
  if we<0.2:
	  we=0.2
#  print n,nr,nr2,nra,nrb,we,float(nrb/50)
  wd=1.
  f2=open(patho+p+'.gp','w')
  f2.write('# Gnuplot script\n\
# crvature/torsion vs arc length of fitted protein axis\n\
reset\n')

  f2.write('set output "'+p+'.ps"\n\
set terminal postscript enhanced portrait 10\n\
set key on\n')
  f2.write('set xtics 10.0\n')
  f2.write('set label font "helvetics,8"\n\
set format x "%4.0f"\n\
set grid\n\n\
### kappa/tau plots ###\n')
  f2.write('set xlabel "res"; unset ylabel\n')
  f2.write('set y2tics mirror\n')

##############################
  f2.write('###1. parameter a plots ###\n\
set ytics 2\n\
set y2tics 0.5\n\
unset y2tics\n\n')

  for i in xrange(nr2):
    b0=100*i
    b1=100*i+50
    b2=100*i+100
    f2.write('set size 1.0,0.475\n\
set origin 0.0,0.5          # parameter a\n')
    if b1>n:
      b1=n
      wd=we
    f2.write('set multiplot\n')
    rsec(p1,b0,b1,f2,wd)
    f2.write('set origin 0.0,0.0          # parameter a\n')
    if b2>n:
      b2=n 
      wd=we
    rsec(p1,b1,b2,f2,wd)
    
    f2.write('unset multiplot\n\n')

  if nra!=0:
    f2.write('set multiplot\n')
    f2.write('set origin 0.0,0.5          # parameter a\n')
    rsec(p1,b2,n,f2,we)
    f2.write('unset multiplot\n\n')
#    f2.write('set multiplot\n')
    
#    f2.write('set xrange ['+str(b0)+':'+str(b1)+']\n'
#    f2.write('set yrange [0:10]\n')
#    f2.write('plot "'+p+'.pa" i 0 u 1:($1 >='+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0) \\\n\
#	  w l lt 1 lw 1 title "{a}", \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0):3 \\\n\
#	  w labels rotate left offset 0,1 notitle, \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0):1 \\\n\
#	  w labels rotate left offset 0,2 notitle, \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b0)+' && $1 <= '+str(b1)+' ? $2 : 1/0) \\\n\
#	  pointtype 7 pointsize 0.5 notitle\n\n')
	  
#    f2.write('set origin 0.0,0.0          # parameter a\n')
#    f2.write('set xrange ['+str(b1)+':'+str(b2)+']\n'
#    f2.write('set yrange [0:10]\n')   
#    f2.write('plot "'+p+'.pa" i 0 u 1:($1 >='+str(b1)+' && $1 <= '+str(b2)+' ? $2 : 1/0) \\\n\
#	  w l lt 1 lw 1 title "{a}", \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b1)+' && $1 <= '+str(b2)+' ? $2 : 1/0):3 \\\n\
#	  w labels rotate left offset 0,1 notitle, \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b1)+' && $1 <= '+str(b2)+' ? $2 : 1/0):1 \\\n\
#	  w labels rotate left offset 0,2 notitle, \\\n\
#	  "'+p+'.pa" i 0 u 1:($1 >= '+str(b1)+' && $1 <= '+str(b2)+' ? $2 : 1/0) \\\n\
#	  pointtype 7 pointsize 0.5 notitle\n\n')

###############################
#2.write spline fit of axis skt
#  print hp
  f2.write('set xlabel "s"; unset ylabel\n\
set xtics 10.0\n\
set ytics 0.5\n\
unset x2tics\n\
set xtics mirror\n')

  f2.write('set label font "helvetics,8"\n\
set format x "%4.0f"\n\
set format y "%4.1f"\n\
set format y2 "%4.1f"\n\
set grid\n\n\
### kappa/tau plots ###\n')
  f2.write('set multiplot\n')
  sal=sa[-1]
  wt=50
#  print len(sa)
  nh=len(hp)
  nh2=nh/2
  nh2a=nh%2
  cog=0
  for h in hp:
 # for i in xrange()
 #   cog+=1
    h0,h1=h[0],h[1]
    a0=h0-1
    a1=h1-1
    b0,b1=a0-1,a1-1
    r0,r1=resi[a0],resi[a1]
#    print h0,h1,r0,r1,a0,a1
    c0,c1=b0*20,b1*20
    s0=sa[c0]
    s1=sa[c1]
#    print s0,s1
    sl=s1-s0
    n=int(sl/wt)
    se=sl-wt*n
#    se=int(round(sl-wt*n))
#    ld=int(round(sl))
    we=float(se)/wt

    for i in xrange(n):
      cog+=1
      la=wt*i+s0
      lb=wt*(i+1)+s0
      if cog%2==1:
        og='0.5'
        gseck(p1,la,lb,f2,1.,og)
      else:
        og='0.0'
        gseck(p1,la,lb,f2,1.,og)
        f2.write('unset multiplot\n\n')
        f2.write('set multiplot\n')


    la=wt*n+s0
    lb=s1
#    print la,lb,we
    if we<0.3:
      we=0.3
#    print la,lb,we
    cog+=1
    if cog%2==1:
      og='0.5'
      gseck(p1,la,lb,f2,we,og)
    else:
      og='0.0'
      gseck(p1,la,lb,f2,we,og)
      f2.write('unset multiplot\n\n')
      f2.write('set multiplot\n')
  f2.write('unset multiplot\n\n')
#3D box

  f2.write('	### 3D axis ###\n\
  set origin 0,0; set size 1,1\n\
  set autoscale\n\
  set format "%2.0f"\n\
  set xtics autofreq; set ytics autofreq; set ztics autofreq\n\
  set xlabel "x"; set ylabel "y"; set zlabel "z"\n')
  f2.write('splot "'+p1+'.axsp" i 0 u 5:6:7 w l lt 1 title "r(s)", \\\n\
        "'+p1+'.axsp" i 0 u 5:6:7:10 w labels notitle\n')

  f2.write('unset multiplot\n\n\
	### 3D backbone & helix axis ###\n\
	set origin 0,0; set size 1,1\n\
	set autoscale\n\
	set format "%2.0f"\n\
	set xtics autofreq; set ytics autofreq; set ztics autofreq\n\
	set xlabel "x"; set ylabel "y"; set zlabel "z"\n\
	set border 4095\n\
splot "'+p+'.dca" i 0 u 3:4:5 w l lt 1 title "r(s)", \\\n\
	      "'+p+'.dca" i 0 u 3:4:5:2 w labels notitle,\\\n\
 "'+p+'.dca" i 0 u 9:10:11 \\\n\
	  pointtype 7 pointsize 0.5 lc rgb "red" notitle\n')

#3D individual helices
  for a in hp:
#    print a
    r3d(p,a[0],a[1],f2)
