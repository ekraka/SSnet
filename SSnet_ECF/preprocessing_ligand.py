import numpy as np
import sys

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

d = get_inps()

f = open(d['data_file'], 'r')
lines = f.readlines()
f.close()

smiles, targets, score = [], [], []

for line in lines:
	if len(line.strip().split()) == 0:
		continue
	a,b,c = line.strip().split()
	smiles.append(a)
	targets.append(b)
	score.append(c)

np.save('ligand_data.npy', {'smiles':smiles, 'targets': targets, 'score': score})
