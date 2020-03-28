import math
import os
import numpy as np
from scipy.signal import find_peaks
import sys
import matplotlib.pyplot as plt

def distance(a,b):
        res=0
        for i in range (len(a)):
                res+=(a[i]-b[i])**2
        return math.sqrt(res)

def job(file, mod = 0):
    if mod == 0:
        f = open(file, 'r')
        lines = f.readlines()
        f.close()
    else:
        lines = file.split('\n')

    b_info = []
    cords = []

    for line in lines:
            k = line.strip().split()
            if len(k) < 1:
                    continue
            if k[0] != 'ATOM':
                    continue
            #print (line)
            x, y, z = list(map(float, [line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]))
            at = line[11:17].strip()
            if at == 'CA':
                    b_info.append(float(k[-1]))
                    cords.append([x,y,z])

    #print (len(b_info))

    peak_ids, _ = find_peaks(b_info)

    temp_lis = []
    for i in peak_ids:
        if b_info[i] > 75.0:
            temp_lis.append(i)

    peak_ids = temp_lis

    clusters = {}
    vs = set()

    for i in peak_ids:
            if i in vs:
                    continue
            clusters[i] = [cords[i]]
            vs.add(i)
            for j in peak_ids:
                    if i == j:
                            continue
                    if distance(cords[i], cords[j]) < 15:
                            vs.add(j)
                            clusters[i].append(cords[j])
    COM = []
    for i in clusters:
            tarr = np.array(clusters[i])
            COM.append([tarr[:,0].mean(), tarr[:,1].mean(), tarr[:,2].mean()])

    return COM

if __name__ == '__main__':
    job(sys.argv[1])

