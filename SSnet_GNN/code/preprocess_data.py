from collections import defaultdict
import os
import pickle
import numpy as np
import wget
import sys
import haxis
from rdkit import Chem


def create_atoms(mol):
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    atoms = [atom_dict[a] for a in atoms]
    return np.array(atoms)


def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict


def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""

    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]

    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict

        for _ in range(radius):

            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                fingerprints.append(fingerprint_dict[fingerprint])
            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    edge = edge_dict[(both_side, edge)]
                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)


def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)


def split_sequence(sequence, ngram):
    sequence = '-' + sequence + '='
    words = [word_dict[sequence[i:i+ngram]]
             for i in range(len(sequence)-ngram+1)]
    return np.array(words)


def dump_dictionary(dictionary, filename):
    with open(filename, 'wb') as f:
        pickle.dump(dict(dictionary), f)

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


if __name__ == "__main__":

    data_file = '../dataset/' + DATASET + '/original/data.txt'
    path = './targets/'


    f = open(d['data_file'], 'r')
    lines = f.readlines()
    f.close()

    os.chdir(path)

    pro_dic = {}

    for line in lines:
            if len(line.strip().split()) == 0:
                    continue
            pdb = line.strip().split()[1] + '.pdb'

            protein = job(pdb)
            if protein is not None:
                    pro_dic[pdb[:-4]] = protein

    DATASET, radius, ngram = sys.argv[1:]
    radius, ngram = map(int, [radius, ngram])

    with open('../dataset/' + DATASET + '/original/data.txt', 'r') as f:
        data_list = f.read().strip().split('\n')

    """Exclude data contains '.' in the SMILES format."""
    data_list = [d for d in data_list if '.' not in d.strip().split()[0]]
    N = len(data_list)

    atom_dict = defaultdict(lambda: len(atom_dict))
    #print (atom_dict)
    #exit()
    bond_dict = defaultdict(lambda: len(bond_dict))
    fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    edge_dict = defaultdict(lambda: len(edge_dict))
    #word_dict = defaultdict(lambda: len(word_dict))

    Smiles, compounds, adjacencies, proteins, interactions = '', [], [], [], []

    for no, data in enumerate(data_list):

        print('/'.join(map(str, [no+1, N])))

        smiles, target, interaction = data.strip().split()
        #print (smiles)

        if target in pro_dic:

            proteins.append(pro_dic[target])
        else:
            print ('Warning:',target, 'is ignored due to some error!')
            continue
        
        Smiles += smiles + '\n'

        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))  # Consider hydrogens.
        atoms = create_atoms(mol)
        i_jbond_dict = create_ijbonddict(mol)

        #print (atom_dict)
        #exit()

        fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius)
        compounds.append(fingerprints)

        adjacency = create_adjacency(mol)
        adjacencies.append(adjacency)

        interactions.append(np.array([float(interaction)]))

    dir_input = ('../dataset/' + DATASET + '/input/'
                 'radius' + str(radius) + '_ngram' + str(ngram) + '/')
    os.makedirs(dir_input, exist_ok=True)

    with open(dir_input + 'Smiles.txt', 'w') as f:
        f.write(Smiles)
    np.save(dir_input + 'compounds', compounds)
    np.save(dir_input + 'adjacencies', adjacencies)
    np.save(dir_input + 'proteins', proteins)
    np.save(dir_input + 'interactions', interactions)
    dump_dictionary(fingerprint_dict, dir_input + 'fingerprint_dict.pickle')
    #dump_dictionary(word_dict, dir_input + 'word_dict.pickle')

    print('The preprocess of ' + DATASET + ' dataset has finished!')
