import wget
import os
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Avalon.pyAvalonTools as A
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

def latent_space(smiles, N_BITS=512):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)
    m = Chem.AddHs(m)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=N_BITS))
    return np.array(A.GetAvalonFP(m))

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

def predic(model, X, proteins, sm, refe_f):
    if not refe_f:
        x = [proteins.reshape((1, 9000, 2)), X.reshape((1, 512))]
    else:
        x = [proteins.reshape((-1, 9000, 2)), X.reshape((-1, 512))]
    ans = model.predict(x)
    l = [[float(ans[i][0]), sm[i]] for i in range (len(sm))]
    l.sort(reverse = True)
    if not refe_f:
        print('SSnet probability:', l[0][0])
    g = open('results_'+sys.argv[1].split('.')[0]+'.txt', 'a')
    for i in range (len(l)):
        g.write(l[i][1]+' '+str(l[i][0])+'\n')
    g.close()
    print ('Done!')

def job(pdb, smile, model_weights):
    if pdb not in os.listdir('.'):
        print ('PDB file not found, trying to download!')
        try :
            url = 'https://files.rcsb.org/download/'+pdb[:4]+'.pdb'
            wget.download(url,pdb[:4]+'.pdb')
        except HTTPError:
            print ('Error while downloading:', pdb[:4])


    filename = pdb.split('.')[0]

    d = mtt.get_kt(pdb, {})

    #np.save('temp.npy', d)
    proteins = np.array(t_k(filename, d)).T
    refe_f = 0
    if smile[-4:] == '.smi':
            refe_f = 1 # reference file 
            X, smile_arr = get_dat_sm(smile)
    else:
            X = latent_space(str(smile))
            smile_arr = [smile]

    Proeins_shape = proteins.shape
    Drags_shape=512#X.shape[0]
    #print (proteins.shape)
    print(Proeins_shape)
    print(Drags_shape)

    multiBranch_Conv_model = create_multiBranch_Conv_model(Proeins_shape,Drags_shape,)
    multiBranch_Conv_model.load_weights(model_weights)
    if refe_f:
        proteins = np.array([proteins for i in range (len(X))])
    #print (proteins_.shape)
    predic(multiBranch_Conv_model, X, proteins, smile_arr, refe_f)
    

if __name__ == '__main__':
    model_weights = os.path.join(os.path.split(sys.argv[0])[0], 'model.h5')
    job(sys.argv[1], sys.argv[2], model_weights)
    


