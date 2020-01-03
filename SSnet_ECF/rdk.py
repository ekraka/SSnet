import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import copy


def latent_space(smiles, N_BITS=512):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)
    m = Chem.AddHs(m)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=N_BITS))

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

d_inp = get_inps()

X, protein, score = [], [], []
count_er = 0

d = np.load('ligand_data.npy').item()
smiles, target, ic50 = d['smiles'], d['targets'], d['score']
for i in range (len(target)):
    if type(target[i]) == list:
        target[i] = str(target[i][0])
    else:
        target[i] = str(target[i].split(',')[0])


t_k_data = np.load(d_inp['targets_dir'] + 'target_data.npy').item()
#print (t_k_data['5tgz'].shape)
def t_k(tar, dat):
    if tar not in dat:
        #print (tar)
        return None
    return dat[tar]

pd = copy.copy(target)
pdbs = []

for i in range (len(target)):
    target[i] = t_k(target[i], t_k_data)


for i in range (len(smiles)):
    #print (smiles[i])
    if 0 and str(smiles[i]) in ['C[C@H]1CN2[C@@H]([C@@H](C)O1)C1(Cc3cc4c(noc4c(F)c23)C2=[N](C)N=NS2)C(=O)NC(=O)NC1=O',
				'CC(C)(C)OC(=O)N(C1CC1)C1=NC(=C[N]2=C(\C=C3/NC(=O)NC3=O)C=NC12)c1cccc(OC(F)(F)F)c1',
				'O=C1NC(=O)\C(N1)=C\C1=[N]2C=C(N=C(NC3CC3)C2N=C1)C#Cc1ccccc1',
				'O=C1NC(=O)\C(N1)=C\C1=[N]2C=C(N=C(NC3CC3)C2N=C1)c1cccc(CN2CCOCC2)c1',
				'FC(F)(F)Oc1cccc(c1)C1=C[N]2=C(\C=C3/NC(=O)NC3=O)C=NC2C(NC2CC2)=N1',
				'Fc1cccc(c1)C1=C[N]2=C(\C=C3/NC(=O)NC3=O)C=NC2C(NC2CC2)=N1',
				'CC(C)(c1nc(nc2N3CCOC[C@H]3COc12)[N]1=C(N)Nc2ccccc12)S(C)(=O)=O',
				'CC(C)(c1nc(nc2N3CCOC[C@H]3COc12)[N]1=C(Nc2ccccc12)N1CCOCC1)S(C)(=O)=O',
				'ONC(=O)CCCCn1cc(C(=O)Nc2cc(nn2)-c2ccccc2)c(=O)c2ccccc12',
				'[O-][N]1=CC(=CC=C1)c1cnc(NCc2c3CCOc3ccc2F)n2cnnc12',
				'c1ccc(c(c1)C=N(=Cc2ccccc2O)N(CCO)C3=[NH+]CCN3)O',
				'c1cc(ccc1C=N(=Cc2ccc(cc2)Cl)N(CCO)C3=[NH+]CCN3)Cl',
				'c1ccc(cc1)C2=CC(=N(=C(C#N)NC(=O)c3ccccc3)C(=C2)c4ccccc4)c5ccccc5',
				'F[Si-2](F)(F)(F)(F)F.[Na+].[Na+]',
				'[NH4+].[NH4+].F[Si-2](F)(F)(F)(F)F',
				'O=Cl=O']:
        continue
    try:
        try:
            x = latent_space(str(smiles[i]))
            if x is None:
                count_er += 1
                continue
        except KeyError:
            count_er += 1
            continue
    except IndexError:
        count_er += 1
        continue

    if target[i] is None:
        #print('protein error',pd[i], ic50[i])
        count_er += 1
        continue

    score.append(ic50[i])
    protein.append(target[i])
    X.append(x)
    pdbs.append(pd[i])

print(count_er)
np.save('final_data.npy', {'X':X, 'proteins':protein, 'score':score, 'pr_names':pdbs})
