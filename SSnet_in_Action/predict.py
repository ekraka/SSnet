import pickle
import os
import numpy as np
import sys
from keras.models import Sequential,Model
import tensorflow as tf
from keras.layers import Dense, Activation, Input,RepeatVector,Embedding, Flatten, Concatenate,Dropout
from keras.models import Model
from keras.utils.vis_utils import model_to_dot
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import average, concatenate,RepeatVector,Lambda,add,subtract
from keras.layers.normalization import BatchNormalization
from keras.utils.vis_utils import model_to_dot
import random
from keras import backend as K
from keras.callbacks import (ModelCheckpoint, LearningRateScheduler,
                             EarlyStopping, ReduceLROnPlateau,CSVLogger)
from keras.layers import Conv2D, MaxPooling2D,Conv1D,GlobalMaxPooling1D,MaxPooling1D,Reshape,Add
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D,AveragePooling1D
import make_target as mtt
from keras.regularizers import l2 
import preprocessing_ligand as ppl 
import grad

# remove if using saved fingerprints
from rdkit import Chem
from rdkit.Chem import AllChem

def conv_blocks(ft_number,k_size,input_tensor):
    x = Conv1D(filters=ft_number, 
                     kernel_size=k_size, data_format='channels_last',
                     padding='same',
                     kernel_regularizer=l2(l2_lambda)
              )(input_tensor)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)

    
#     x = Conv1D(filters=ft_number, 
#                      kernel_size=k_size, 
#                      padding='same',
#                      kernel_regularizer=l2(l2_lambda)
#               )(x)
#     x = Activation('relu')(x)
#     x = BatchNormalization()(x)
    
    return x
def op(inputs):
    x, y = inputs
    return K.pow((x - y), 2)

def conv_branch(init_input,kernel_size):
    x=conv_blocks(ft_number=32,k_size=kernel_size,input_tensor=init_input)
    #x=MaxPooling1D(2,padding='same')(x)

    x=conv_blocks(ft_number=64,k_size=kernel_size,input_tensor=x)
    #x=MaxPooling1D(2,padding='same')(x)

    x=conv_blocks(ft_number=128,k_size=kernel_size,input_tensor=x)
    u = GlobalMaxPooling1D()(x)
    u_broadcast=RepeatVector(x.shape[1])(u)

    o=Lambda(op)([u_broadcast,x])  # K.pow((x - y), 2) 
    var = GlobalMaxPooling1D()(o)
    X_vector = concatenate([u,var])
    
    #X_vector = Dense(64)(X_vector)
    return X_vector

l2_lambda=0.05
def create_multiBranch_Conv_model(Proeins_shape,Drags_shape):
    
    # left branch --> dimension reduction for Proeins_shape
    
    proeins_input_tensor = Input(shape=Proeins_shape, name='proeins_input_tensor')
    #init_input = Reshape((Proeins_shape[1], Proeins_shape[0]),input_shape=Proeins_shape,name='init_input')(proeins_input_tensor)
    init_input = proeins_input_tensor
    #branch 0
    w=conv_branch(init_input,10)
    #branch 1
    x=conv_branch(init_input,5)

    #branch 2
    y=conv_branch(init_input,15)
    
    #branch 3
    z=conv_branch(init_input,20)
    
    protein_concat = concatenate([w,x,y,z], name='protein_concat_')
    
    protein_dense = Dense(128)(protein_concat)
    protein_dense = BatchNormalization()(protein_dense)
    protein_dense = Activation('relu')(protein_dense)
    protein_dense = Dropout(0.5)(protein_dense)
    
    
    # right branch --> dimension reduction for drug/ligand
    drag_input_tensor = Input(shape=(Drags_shape,),name='drag_input_tensor')
    d = Dense(128)(drag_input_tensor)
    d = BatchNormalization()(d)
    d = Activation('relu')(d)
    d = Dropout(0.5)(d)

    
#     # merge the branches together
    final_branch = concatenate([protein_dense,d], name='protein_darg_concat_')
    
    final_dense = Dense(64)(final_branch)
    final_dense = BatchNormalization()(final_dense)
    final_dense = Activation('relu')(final_dense)
    final_dense = Dropout(0.5)(final_dense)

#     final_dense = Dense(32)(final_dense)
#     final_dense = BatchNormalization()(final_dense)
#     final_dense = Activation('relu')(final_dense)
#     final_dense = Dropout(0.5)(final_dense)
    
    final_output = Dense(1, activation='sigmoid', name='final_output')(final_dense)
    
#     model = Model(inputs=proeins_input_tensor, outputs=protein_concat)
    model = Model(inputs=[proeins_input_tensor,drag_input_tensor], outputs=final_output)
    return model

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

def latent_space(smiles, N_BITS=512):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)
    m = Chem.AddHs(m)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=N_BITS))
    #return np.array(A.GetAvalonFP(m))

def round_sig(x, sig=3):
    if x < 1e-2:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)


def get_dat_sm(smiles):
    f= open(smiles, 'r')
    lines = f.readlines()
    f.close()

    sm, X = [], []
    sm_dic = {}
    count_er = 0 
    for line in lines:
        k = line.strip().split()
        if len(k) > 0:
            s = k[0]
            if s not in sm_dic:
                x = latent_space(s)
                sm_dic[s] = x
            else:
                x = sm_dic[s]
            if x is None:
                count_er += 1
                continue
            sm.append(s)
            X.append(x)
    print (count_er)
    return np.array(X), sm

def predic(model, X, proteins, sm, refe_f, tar = None, sm_name = None):
    global parms
    if not refe_f:
        x = [proteins.reshape((1, 9000, 2)), X.reshape((1, 512))]
    else:
        x = [proteins.reshape((-1, 9000, 2)), X.reshape((-1, 512))]
    ans = model.predict(x)
    if tar is not None:
        l = [[float(ans[i][0]), sm[i], tar[i]] for i in range (len(sm))]
    else:
        l = [[float(ans[i][0]), sm[i]] for i in range (len(sm))]
    
    l.sort(reverse = True)
    if not refe_f:
        print('SSnet probability:', l[0][0])
    if sm_name:
        g = open('results_'+parms['-t'].split('.')[0]+'_' + parms['-l'].split('.')[0] + '.txt', 'w')
    elif '-i' in parms:
        g = open('results_'+parms['-i'].split('.')[0] +'.txt', 'w')
    else:
        g = open('results_'+parms['-t'].split('.')[0] +'.txt', 'w')
    for i in range (len(l)):
        if tar is not None:
            g.write(l[i][1]+' '+str(l[i][2])+' '+str(l[i][0])+'\n')
        else:
            g.write(l[i][1]+' '+str(l[i][0])+'\n')
    g.close()
    print ('Done!')

def get_pdb_data(pdb, targets_dir = '.'):
    if pdb not in os.listdir(targets_dir):
        import wget
        print ('PDB file not found, trying to download!', pdb[:-4])
        try :
            url = 'https://files.rcsb.org/download/'+pdb[:-4]+'.pdb'
            wget.download(url,os.path.join(targets_dir,pdb))
        except:
            print ('Error while downloading:', pdb[:-4])
            return None


    filename = pdb[:-4]

    d = mtt.get_kt(os.path.join(targets_dir, pdb), {})

    #np.save('temp.npy', d)
    return np.array(t_k(filename, d)).T

def get_sm_tar_data(smiles, targets, targets_dir):
    X, proteins, P = [], [], []
    sm , tar = [], []
    sm_dic = {}
    count_er = 0

    tar_dic = {}

    for i in range (len(smiles)):
        s = smiles[i]
        if s not in sm_dic:
            x = latent_space(s)
            sm_dic[s] = x
        else:
            x = sm_dic[s]
        if x is None:
            count_er += 1
            continue
        try:
            #print (os.path.join(targets_dir,targets[i] + '.pdb'))
            if targets[i] not in tar_dic:
                p = get_pdb_data(targets[i] + '.pdb', targets_dir)
                #print (p)
                tar_dic[targets[i]] = p 
            else:
                p = tar_dic[targets[i]]
        except TypeError:
            count_er += 1
            continue 
        sm.append(s)
        tar.append(targets[i])
        proteins.append(p)
        X.append(x)
    print ('Number of errors', count_er)
    return np.array(sm), np.array(X), np.array(tar), np.array(proteins)


def pro_from_file(input_file, targets_dir):
    smiles, targets = ppl.job(input_file) # provide single input as inp
    #print (smiles, targets)
    sm, X, tar, proteins = get_sm_tar_data(smiles, targets, targets_dir)
    return sm, X, tar, proteins


    


def job(pdb, smile, model_weights, file_ref = 1):
    #print (file_ref)
    if file_ref < 4 or file_ref > 5:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    Proeins_shape = (9000,2)
    Drags_shape=512#X.shape[0]
    multiBranch_Conv_model = create_multiBranch_Conv_model(Proeins_shape,Drags_shape,)
    multiBranch_Conv_model.load_weights(model_weights)

    if file_ref == 1: # file for PLI
        sm, X, tar, proteins = pro_from_file(pdb, smile) # input_file, proteins_directory
        #print (proteins)
        predic(multiBranch_Conv_model, X, proteins, sm, 1, tar)

    elif file_ref == 2: # for a pdb and .smi 
        X, smile_arr = get_dat_sm(smile)
        proteins = get_pdb_data(pdb)
        proteins = np.array([proteins for i in range (len(X))])
        predic(multiBranch_Conv_model, X, proteins, smile_arr, 1, None, 1)
    elif file_ref == 3: # a protein and a smiles
        X = latent_space(str(smile))
        smile_arr = [smile]
        proteins = get_pdb_data(pdb)
        predic(multiBranch_Conv_model, X, proteins, smile_arr, 0)
    elif file_ref == 4: # from saved fingerprints as .npy
        data = np.load(smile).item()
        X = np.array(data['X'])
        proteins = get_pdb_data(pdb)
        proteins = np.array([proteins for i in range (len(X))])
        predic(multiBranch_Conv_model, X, proteins, data['smile_arr'], 1, None, 1)
    elif file_ref == 5: # from saved fingerprints as .dat
        with open(smile, 'rb') as infile:
            data = pickle.load(infile)
        X = np.array(data['X'])
        proteins = get_pdb_data(pdb)
        proteins = np.array([proteins for i in range (len(X))])
        predic(multiBranch_Conv_model, X, proteins, data['smile_arr'], 1, None, 1)
    elif file_ref == 6: # protein and smile for grad
        X = latent_space(str(smile))
        smile_arr = [smile]
        proteins = get_pdb_data(pdb)
        grad.job(proteins, X, pdb, multiBranch_Conv_model)
        #predic(multiBranch_Conv_model, X, proteins, smile_arr, 0)

def pars_parm():
    print ("""
======================================================================
==       Secondary Structure Based Deep Neural Network Model        ==
==          for Protein-ligand Interaction Prediction :             ==
==                         Code version 1.0                         ==
==                                                                  ==
==    Computational and Theoretical Chemistry Group (CATCO), SMU    ==
==                     Dallas, Texas 75275 USA                      ==
======================================================================


            _|_|_|    _|_|_|                        _|      
            _|        _|        _|_|_|    _|_|    _|_|_|_|  
            _|_|      _|_|    _|    _|  _|_|_|_|    _|      
                _|        _|  _|    _|  _|          _|      
            _|_|_|    _|_|_|  _|    _|    _|_|_|      _|_|  
                                                            
If no idea how it works!, include -h as argument                                                           
                                                            
                                                            """)
    # Default
    d = {'-m': '10', '-t_dir': '.'}
    for i in range (1, len(sys.argv)):
        if sys.argv[i] == '-h':
            d['-h'] = 1
        elif '-' in sys.argv[i]:
            d[sys.argv[i]] = sys.argv[i+1]
    return d 

def help_args():
    print("""
Check Github page for details:
https://github.com/ekraka/SSnet 

SSnet can be used with one of the four different models:
-m <mode>
    i) mode = 10 (Default)
        (This model was trained based on 10nM IC50 cutoff for actives)
    ii) mode = 25
        (This model was trained based on 25nM IC50 cutoff for actives)
    iii) mode = 100
        (This model was trained based on 100nM IC50 cutoff for actives)
    iv) mode = grad
        (This model was trained for generating heatmaps for potential binding location)

Multiple proteins and ligands can be provided as input file
-i <file>

Single protein can be described as
-p <pdb_file>

ligands can be provided as either SMILES or a .smi file with multiple ligands
-l <smi>
    
    """)
    

if __name__ == '__main__':

    parms = pars_parm()
    print ('Input parameters:', parms)
    
    if '-m' in sys.argv:
        m = parms['-m']
        if m == '10':
            model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'models/model_10.h5')
        elif m == '100':
            model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'models/model_100.h5')
        elif m == '25':
            model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'models/model_25.h5')
        elif m == 'grad':
            model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'models/model_grad.h5')
        else:
            raise Exception('-m mode unknown!')
    else:
        model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'models/model_10.h5')

    if '-i' in parms:
        job(parms['-i'], parms['-t_dir'], model_weights, 1)
    elif '-h' in parms:
        help_args()
        print ('Help might come through others ways too! Email: nirajverma288@gmail.com')
    else:
        target, ligand = parms['-t'], parms['-l']
        if ligand[-4:]  == '.smi':
            job(parms['-t'], parms['-l'], model_weights, 2)
        elif ligand[-4:]  == '.npy':
            job(parms['-t'], parms['-l'], model_weights, 4)
        elif ligand[-4:]  == '.dat':
            job(parms['-t'], parms['-l'], model_weights, 5)
        elif parms['-m'] == 'grad':
            job(parms['-t'], parms['-l'], model_weights, 6)
        else:
            job(parms['-t'], parms['-l'], model_weights, 3)
    



