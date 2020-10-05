import os
import haxis
import numpy as np
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Avalon.pyAvalonTools as A
import make_box as m_box
import pandas as pd
from keras.models import Sequential,Model
import tensorflow as tf
from keras.layers import Dense, Activation, Input,RepeatVector,Embedding, Flatten, Concatenate,Dropout
from keras.models import Model
from keras.utils.vis_utils import model_to_dot
from sklearn import metrics as mt
from keras import metrics
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import average, concatenate,RepeatVector,Lambda,add,subtract
from keras.layers.normalization import BatchNormalization
from keras.utils.vis_utils import model_to_dot
from sklearn.utils import class_weight, shuffle
import random
from keras import backend as K
from keras.callbacks import (ModelCheckpoint, LearningRateScheduler,
                             EarlyStopping, ReduceLROnPlateau,CSVLogger)
from sklearn.metrics import precision_recall_fscore_support, classification_report
import pickle
from keras.layers import Conv2D, MaxPooling2D,Conv1D,GlobalMaxPooling1D,MaxPooling1D,Reshape,Add
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D,AveragePooling1D
from keras.regularizers import l2 
from sklearn.metrics import roc_auc_score

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

def get_heatm(conv, x, model, ref = 0):
    #model = multiBranch_Conv_model

    #print (model.layers)

    #x = [train_proteins[index].reshape((-1, 2, 9000)), train_X[index].reshape((-1, 512))]

    predicted_class = int(round(model.predict(x)[0][0]))

    #predicted_class_output = model.get_layer('protein_concat_').output #model.output[0]
    predicted_class_output = model.output#[:, predicted_class]
    
    #print (predicted_class_output)
    if ref:

            pass

    last_conv_layer = model.get_layer(conv)

    # This is the gradient of the predicted class with regard to
    # the output feature map of `block5_conv3`
    grads = K.gradients(predicted_class_output, last_conv_layer.output)[0]
    
    #print (grads.shape)
    a, s = grads.shape[-2:]
    #if type(a) != int:
    #    a = 1
    #print (a, s)

    # This is a vector of shape (512,), where each entry
    # is the mean intensity of the gradient over a specific feature map channel
    ##pooled_grads = grads
    pooled_grads = K.mean(grads, axis=(0,1))
    
    #print (pooled_grads.shape)

    # This function allows us to access the values of the quantities we just defined:
    # `pooled_grads` and the output feature map of `block5_conv3`,
    # given a sample image
    #print (last_conv_layer.output.shape)
    iterate = K.function([model.input[0], model.input[1]], [pooled_grads, last_conv_layer.output[0]])

    # These are the values of these two quantities, as Numpy arrays,
    # given our sample image

    #print (iterate(x))
    pooled_grads_value, conv_layer_output_value = iterate(x)
    
    #print (conv_layer_output_value.shape, pooled_grads_value.shape)

    #pooled_grads_value = pooled_grads_value.reshape((conv_layer_output_value.shape))
    
    


    # We multiply each channel in the feature map array
    # by "how important this channel is" with regard to the predicted class
    for i in range(s):
        conv_layer_output_value[:,i] *= pooled_grads_value[i]

    # The channel-wise mean of the resulting feature map
    # is our heatmap of class activation
    #print (conv_layer_output_value.shape)
    heatmap = np.mean(conv_layer_output_value, axis=(-1,))
    heatmap = np.maximum(heatmap, 0)
    if np.max(heatmap) > 0:
        heatmap /= np.max(heatmap)
    return heatmap

def shape_heatmap(heatmap):
    f = int(9000 / heatmap.shape[0])
    print ('Factor:', f)
    l = []
    for i in heatmap:
        l += [i]*f
    if len(l) < 9000:
        l += [0]*(9000 - len(l))
    return np.array(l).reshape((6, 1500))

from math import log10, floor
def round_sig(x, sig=3):
    if x < 1e-2:
        return 0.0
    return round(x, sig-int(floor(log10(abs(x))))-1)

def create_pdb(pdb, heatmap):
    gami, gama = np.amin(heatmap), np.amax(heatmap)
    for i in range (len(heatmap)):
        a = heatmap[i]
        #gami, gama = a.min(), a.max()
        heatmap[i] = np.interp(a, (gami, gama), (0.0, 100.0))
        
        #print (np.max(heatmap))
        
    f = open(pdb, 'r')
    lines = f.readlines()
    f.close()
    
    st, ref = [], []
    mod = 0 # model
    res = [] # residue
    for line in lines:
        if line.strip().split()[0] in ['MODEL', 'ENDMDL']:
            st.append(line.strip())
            ref.append(None)
        if 'TER' in line[:3]:
            st.append(line.strip())
            ref.append(None)
            mod += 1
            res = []
        if 'ATOM' in line[:4]:
            nres = int(line[22:26])
            if nres not in res:
                res.append(nres)
            st.append(line[:61])
            ref.append(round_sig(heatmap[mod][len(res) - 1]))
    print ('Number of chains:', mod)
    k_ref = []
    for i in ref:
        if i is not None:
            k_ref.append(i)
        else:
            k_ref.append(ref[0])
            
    k_ref = np.array(k_ref)
    k_ref = np.interp(k_ref, (k_ref.min(), k_ref.max()), (0.0, 100.0))
    
    fst = ''
    
    for i in range (len(st)):
        if ref[i] is not None:
            ans = str(round_sig(k_ref[i]))
            #print (ans)
            
            fst += st[i][:-1]+' '*(6-len(ans))+ans + '\n' 
            #print (fst)
            #stop()
        else:
            fst += st[i]+ '\n'
        #print (ref[i])
    '''
    clusters = m_box.job(fst, 1)
    count = i
    for i in clusters:
        print (i)
        l = ['HETATM', str(count), 'Zn', 'Zn', 'A', str(count), round(i[0]), round(i[1]), round(i[2]), '1.00', '0.00', 'Zn']
        fst += "{:>6}{:>5}{:>4} {:>4}{:>2}{:>4}{:>12}{:>8}{:>8}{:>6}{:>6}{:>12}".format(*l)+'\n'
        count += 1
    '''
    g = open('test_'+pdb[-8:], 'w')
    #print (fst)
    g.write(fst)
    g.close()

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

def predic(model, X, proteins, sm):
    x = [proteins.reshape((-1, 9000, 2)), X.reshape((-1, 512))]
    ans = model.predict(x)
    l = [[float(ans[i][0]), sm[i]] for i in range (len(sm))]
    l.sort(reverse = True)
    g = open('results_'+sys.argv[1].split('.')[0]+'.txt', 'w')
    for i in range (len(l)):
        g.write(l[i][1]+' '+str(l[i][0])+'\n')
    g.close()
    print ('Done!')

def job(pdb, smile):
    if pdb not in os.listdir('.'):
        try :
            import wget
            url = 'https://files.rcsb.org/download/'+pdb[:4]+'.pdb'
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
    #np.save('temp.npy', d)
    proteins = np.array(t_k(filename, d)).T
    refe_f = 0
    if smile[-4:] == '.smi':
            refe_f = 1 # reference file 
            X, smile_arr = get_dat_sm(smile)
    else:
            X = latent_space(str(smile))

    Proeins_shape = proteins.shape
    Drags_shape=512#X.shape[0]
    print(Proeins_shape)
    print(Drags_shape)

    #proteins = proteins.reshape((1, 9000, 2))
    #X = X.reshape((1, -1))

    multiBranch_Conv_model = create_multiBranch_Conv_model(Proeins_shape,Drags_shape,)
    multiBranch_Conv_model.load_weights(weight_file)

    #print (multiBranch_Conv_model.summary())

    if refe_f:
            proteins_ = np.array([proteins for i in range (len(X))])
            #print (X.shape, proteins_.shape)
            predic(multiBranch_Conv_model, X, proteins_, smile_arr)
            return
    
    l = ['conv1d_3', 'conv1d_6', 'conv1d_9', 'conv1d_12']
    heatmap = None
    x = [proteins.reshape((-1, 9000, 2)), X.reshape((-1, 512))]
    for i in l:
            if heatmap is None:
                    heatmap = get_heatm(i, x, multiBranch_Conv_model, 1)/len(l)
            else:
                    heatmap += get_heatm(i, x, multiBranch_Conv_model)/len(l)

    heatmap = shape_heatmap(heatmap)
    create_pdb(pdb, heatmap)

if __name__ == '__main__':
    weight_file = '/users/nirajv/scratch/test/SSnet/GRAD/multiBranch_Conv_model_copy.h5'
    job(sys.argv[1], sys.argv[2])
    


