import numpy as np
import wget
import sys
import haxis
import os
from urllib.error import HTTPError

def get_inps():

	f = open(sys.argv[1], 'r')
	lines = f.readlines()
	f.close()

	d = {}
	for line in lines:
		if len(line.strip().split()) > 1:
			a,b = line.strip().split()
			d[a] = b

	if d['script_dir'] == '.':
		d['script_dir'] = ''

	return d

def check_chain(fi):
        f = open(fi, 'r')
        lines = f.readlines()
        f.close()

        d = {}
        for line in lines:
                if 'ENDMDL' in line:
                        break
                if 'ATOM' == line.strip().split()[0]:
                        temp = line[20:22].strip()
                        if temp not in d:
                                d[temp] = 1
        return d


def make(i):
        if i[-4:] != '.pdb':
                return None
        k = check_chain(i)
        if len(k) > 8:
                print (i,'is ignored due to large number of chains !')
                return None

        g = open('temp', 'w')
        g.write(""" 



""")
        name = i.split('.')[0]
        for j in k:
                g.write(name+'\t'+j+'\n')

        g.close()

        return 1

def data(fi):
        f = open(fi, 'r')
        lines = f.readlines()
        f.close()

        k, t = [], []
        ref = 1

        for line in lines:
                if ref and len(line.strip().split()) == 0:
                        break
                if ref:
                        temp = line.strip().split()
                        k.append(float(temp[1]))
                        t.append(float(temp[2]))
        #print fi
        #print lines[12]
        return k,t

def t_k(tar, dat):
    if tar not in dat:
        #print (tar)
        return None
    #print dat[tar][0]
    t, k = [], []
    if len(dat[tar]) > 6:
        #print (tar, len(dat[tar]))
        return None
    for i,j in dat[tar]:
        if len(i) > 1500 or len(i) < 1:
            print (tar, len(dat[tar]), len(i))

            return None
        t += i + [0]*(1500 -len(i))
        k += j + [0]*(1500 -len(i))
        #print len(i)
        #exit()
    return [t+[0]*(9000 - len(t)), k+[0]*(9000 - len(k))]#k+[0]*(18000 - len(t+k)) 

def job(pdb):
    if pdb not in os.listdir('.'):
        try :
            url = 'https://files.rcsb.org/download/'+pdb[:4].upper()+'.pdb'
            wget.download(url,pdb[:4]+'.pdb')
        except HTTPError:
            print ('Error while downloading:', pdb[:4])
            return None

    make(pdb)
    haxis.haxis('temp')
    d = {}

    filename = pdb.split('.')[0]

    for file in os.listdir('.'):
            if file[-5:] == '.htxt' and file[:len(filename)] == filename:
                    k, t = data(file)
                    if file[:-6] in d:
                            d[file[:-6]].append([k,t])
                    else:
                            d[file[:-6]] = [[k,t]]

    proteins = np.array(t_k(filename, d)).T

    return proteins


if __name__ == '__main__':
        d = get_inps()
        data_file = d['data_file']
        path = d['targets_dir']


        f = open(d['data_file'], 'r')
        lines = f.readlines()
        f.close()

        os.chdir(path)

        targets = {}

        for line in lines:
                if len(line.strip().split()) == 0:
                        continue
                pdb = line.strip().split()[1] + '.pdb'

                if pdb[:-4] in targets:
                        continue
                try:
                        protein = job(pdb)
                except ZeroDivisionError:
                        print (pdb, 'ignored due to math error!')
                if protein is not None:
                        targets[pdb[:-4]] = protein

        np.save('target_data.npy', targets)




