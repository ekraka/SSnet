# This script is designed to do a new spline interpolation
# in order to avoid the discontinuity in the torsion plot...  
# 
# Y.T. 2020/02/17
#

import os
import numpy as np
from numpy import linalg as LA
from scipy import interpolate,integrate
from math import sqrt
import sys


def get_alpha_c_xyz(inp_f):

    x=[[]]
    y=[[]]
    z=[[]]
    with open(inp_f) as f1:
      for line in f1:
         if 'ENDMDL' == line.strip().split()[0]:
            break
         if 'TER' == line.strip().split()[0]:
            x.append([])
            y.append([])
            z.append([])
         if 'ATOM' != line.strip().split()[0]:
            continue
         if len(line) > 3:
            id,at,rt,_,_0 = line[7:11].strip(), line[11:17].strip(), line[17:20].strip(), line[20:22].strip(), line[22:27].strip()
            x1,y1,z1 = line[27:38].strip(), line[39:46].strip(), line[47:54].strip()
            if at == 'CA':
                x[-1].append( float(x1) )
                y[-1].append( float(y1) )
                z[-1].append( float(z1) )
    return x,y,z



def find_pairs(a,A):

    pointer = 0
    max_p = len(A[0])

    res_idx = [] # result index

    for i in range(len(a[0])):
       while_true = True
       min_dis = 999.0
       a1 = np.array([a[0][i],a[1][i],a[2][i]])
       while (while_true):
          if (pointer+1) > max_p:
             while_true = False
             break

          A1 =  np.array([A[0][pointer],A[1][pointer],A[2][pointer]])
          dis = LA.norm(a1-A1)
          if dis<= min_dis:

             min_dis = dis
             pointer = pointer + 1
             continue
          else:
              
             res_idx.append( pointer-1 )
             while_true = False
             



    
    res_idx.append( max_p - 1 )

    return res_idx 

def get_init_kt(x,y,z):
   pts_per_residue = 20
   residue = len(x)
   tck, u = interpolate.splprep([x,y,z], k=4, s=0.0)
   x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)


   u_fine = np.linspace(0,1,pts_per_residue*(residue-1))
   x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)


   x_1st, y_1st, z_1st = interpolate.splev(u_fine, tck, der=1, ext=1)

   func = []
   for i in range(len(x_1st)):
       func.append(  sqrt(x_1st[i]*x_1st[i]+y_1st[i]*y_1st[i]+z_1st[i]*z_1st[i]) )
   arc_lengths = integrate.cumtrapz(func,u_fine,initial=0.0) # total arc length...      

   x_2nd, y_2nd, z_2nd = interpolate.splev(u_fine, tck, der=2, ext=1)
   curvatures = []
   for i in range(len(x_2nd)):
       curvatures.append( sqrt( (z_2nd[i]*y_1st[i]-y_2nd[i]*z_1st[i])**2 + \
                                (x_2nd[i]*z_1st[i]-z_2nd[i]*x_1st[i])**2 + \
                                (y_2nd[i]*x_1st[i]-x_2nd[i]*z_1st[i])**2   ) /\
                          (x_1st[i]*x_1st[i]+y_1st[i]*y_1st[i]+z_1st[i]*z_1st[i])**1.5)

   x_3rd, y_3rd, z_3rd = interpolate.splev(u_fine, tck, der=3, ext=1)
   torsions = []
   for i in range(len(x_3rd)):
       torsions.append(  (x_3rd[i]*(y_1st[i]*z_2nd[i]-y_2nd[i]*z_1st[i]) + \
                          y_3rd[i]*(z_1st[i]*x_2nd[i]-z_2nd[i]*x_1st[i]) + \
                          z_3rd[i]*(x_1st[i]*y_2nd[i]-x_2nd[i]*y_1st[i])   ) / \
                        ((y_1st[i]*z_2nd[i]-y_2nd[i]*z_1st[i])**2 + \
                         (z_1st[i]*x_2nd[i]-z_2nd[i]*x_1st[i])**2 + \
                         (x_1st[i]*y_2nd[i]-x_2nd[i]*y_1st[i])**2 ) \
                      )

   alpha_c_idx = find_pairs([x,y,z],[x_fine, y_fine, z_fine])
   k,t = [], []

   #print " "
   #print "   x        y        z        s        k        t        kCA      tCA"
   for i in range(len(x_fine)):
       if i not in alpha_c_idx:
          pass
          #print("% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f" % (x_fine[i], y_fine[i], z_fine[i],arc_lengths[i],curvatures[i],torsions[i]) )
       else:
          #print("% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f% 9.3f" % (x_fine[i], y_fine[i], z_fine[i],arc_lengths[i],curvatures[i],torsions[i],                                                                  curvatures[i],torsions[i]))
          k.append(curvatures[i])
          t.append(torsions[i])
   return k,t


def get_kt(i, d):
   #residue = int(sys.argv[2])#76
   inp_f = i#sys.argv[1]#"c_alpha.txt"

   x,y,z = get_alpha_c_xyz(inp_f)
   d[i[:-4]] = []
   for j in range (len(x)):
       x1, y1, z1 = x[j], y[j], z[j]
       if len(x1) < 5:
           continue
       k, t = get_init_kt(x1, y1, z1)
       d[i[:-4]].append([k,t])

   return d



def job():
    d = {}
    for i in os.listdir('.'):
        if i[-4:] == '.pdb':
            d = get_kt(i, d)
    np.save('target_data.npy', d)
            


if __name__== "__main__":
   job()
   exit()
   








